#include <getopt.h>
#include <sqlite3.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "mbtiles.hpp"
#include "mvt.hpp"
#include "protozero/pbf_reader.hpp"

using namespace std;

bool verbose = false;

struct zxy {
  long long z;
  long long x;
  long long y;
  zxy(long long _z, long long _x, long long _y) {
    z = _z;
    x = _x;
    y = _y;
  }
  bool operator<(zxy const &other) const {
    if (z < other.z) {
       return true;
    }
    if (z > other.z) {
      return false;
    }
    if (x < other.x) {
      return true;
    }
    if (x > other.x) {
      return false;
    }
    if (y < other.y) {
      return true;
    }
    return false;
  }

  bool operator==(zxy const &other) const {
    return z == other.z && x == other.x && y == other.y;
  }

};

void create_tile_set(const string fname, set<zxy>* tile_set) {
  sqlite3 *db;
  if (sqlite3_open(fname.c_str(), &db) != SQLITE_OK) {
    cerr << fname << ": " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  const char *sql = "SELECT zoom_level, tile_column, tile_row from tiles "
      "order by zoom_level, tile_column, tile_row;";

  sqlite3_stmt *stmt;
  if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) != SQLITE_OK) {
    cerr << fname << ": select failed: " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  zxy tile(0, 0, 0);
  while (sqlite3_step(stmt) == SQLITE_ROW) {
    tile.z = sqlite3_column_int(stmt, 0);
    tile.x = sqlite3_column_int(stmt, 1);
    tile.y = sqlite3_column_int(stmt, 2);
    tile_set->insert(tile);
    if (verbose) {
      cout << fname << " " << tile.z << " " << tile.x << " " <<  tile.y
           << " (flipped y: " << (1LL << tile.z) - 1 - tile.y << ")" << endl;
    }
  }

  sqlite3_finalize(stmt);
  if (sqlite3_close(db) != SQLITE_OK) {
    cerr << fname << ": " << "could not close database: "
         << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
}

void classify_tiles(
    const string region_fname, const string boundary_fname,
    const string main_fname, set<zxy>* interior_set,
    set<zxy>* boundary_set, set<zxy>* main_set) {
  set<zxy> region_set;
  create_tile_set(region_fname, &region_set);
  create_tile_set(boundary_fname, boundary_set);
  set_difference(region_set.begin(), region_set.end(),
                      boundary_set->begin(), boundary_set->end(),
                      inserter(*interior_set, interior_set->end()));

  create_tile_set(main_fname, main_set);
  set<zxy> exterior_set;
  set_difference(main_set->begin(), main_set->end(),
                      region_set.begin(), region_set.end(),
                      inserter(exterior_set, exterior_set.end()));

  cout << "Total number of tiles in main tileset: " << main_set->size()
       << endl;
  cout << "Number of interior tiles: " << interior_set->size() << endl;
  cout << "Number of boundary tiles: " << boundary_set->size() << endl;
  cout << "Number of exterior tiles: " << exterior_set.size() << endl;
  cout << "Number of not classified tiles: "
       << main_set->size() - interior_set->size() - boundary_set->size() -
             exterior_set.size() << endl;
}

void delete_interior_tiles(string fname, set<zxy> main_tiles,
                           set<zxy> interior_tiles) {
  sqlite3 *db;
  if (sqlite3_open(fname.c_str(), &db) != SQLITE_OK) {
    cerr << fname << ": " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  const int kNumInteriorTiles = interior_tiles.size();
  int num_tiles_deleted = 0;
  for (auto t : interior_tiles) {
    char query[200];
    snprintf(
        query, 200,
        "delete from tiles where zoom_level=%lld and tile_column=%lld "
        "and tile_row=%lld;",
        t.z, t.x, t.y);
    char *err = 0;
    int rc = sqlite3_exec(db, query, NULL, NULL, &err);
    if (rc != SQLITE_OK){
      cerr << "SQL error: " << err << endl;
      sqlite3_free(err);
    } else {
      num_tiles_deleted++;
      if (num_tiles_deleted == kNumInteriorTiles) {
        cout << "All interior tiles deleted: "
             << num_tiles_deleted << "/" << kNumInteriorTiles << endl;
      } else if(num_tiles_deleted%1000 == 0) {
        cout << "Interior tiles deleted: "
             << num_tiles_deleted << "/" << kNumInteriorTiles << endl;
      }
      if (verbose) {
        cout << "Row " << t.z << "/" << t.x << "/" << t.y
             << " deleted successfully." << endl;
      }
    }
  }

  if (sqlite3_close(db) != SQLITE_OK) {
    cerr << fname << ": could not close database: " << sqlite3_errmsg(db)
         << endl;
    exit(EXIT_FAILURE);
  }
  set<zxy> main_set_after_deletion;
  create_tile_set(fname, &main_set_after_deletion);
  cout << "Number of main tiles after deletion: "
       << main_set_after_deletion.size() << endl;

  // Verify result after deletion.
  set<zxy> main_deleted_set;
  set_difference(main_tiles.begin(), main_tiles.end(),
                      main_set_after_deletion.begin(),
                      main_set_after_deletion.end(),
                      inserter(main_deleted_set, main_deleted_set.end()));
  set<zxy> diff_expected_real;
  set_difference(interior_tiles.begin(), interior_tiles.end(),
      main_deleted_set.begin(), main_deleted_set.end(),
      inserter(diff_expected_real, diff_expected_real.end()));
  if (!diff_expected_real.empty()) {
    cout << "Total number of failed deletions: "
         << diff_expected_real.size() << endl;
    for (auto t : diff_expected_real) {
      cout << "Failed to delete tile  " << t.z << " " << t.x << " "
           << t.y << " (flipped y: " << (1LL << t.z) - 1 - t.y << ")"
           << endl;
    }
  }
}


void read_tile(sqlite3* db, const zxy& t, mvt_tile* tile) {
  const string query = "SELECT tile_data from tiles where zoom_level=? "
      "and tile_column=? and tile_row=?;";
  sqlite3_stmt *stmt;
  if (sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
    cerr << "select failed: " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  sqlite3_bind_int(stmt, 1, t.z);
  sqlite3_bind_int(stmt, 2, t.x);
  sqlite3_bind_int(stmt, 3, t.y);

  string message;
  if (sqlite3_step(stmt) == SQLITE_ROW) {
    int len = sqlite3_column_bytes(stmt, 0);
    const char *s = (const char *) sqlite3_column_blob(stmt, 0);
    message = string(s, len);
  }
  sqlite3_finalize(stmt);

  bool was_compressed;
  try {
    if (!tile->decode(message, was_compressed)) {
      cerr << "Couldn't parse tile " << t.z << " " << t.x << " "
           << t.y << endl;
      exit(EXIT_FAILURE);
    }
  } catch (protozero::unknown_pbf_wire_type_exception e) {
    cerr << "PBF decoding error in tile " << t.z << " " << t.x << " "
         << t.y << endl;
    exit(EXIT_FAILURE);
  }
}

// TODO: Implement this function.
bool feature_intersects_region(const mvt_feature& region_feature,
                               const mvt_feature& main_feature) {
  return true;
}

int clear_names(mvt_layer* layer, mvt_feature* feature) {
  int num_names_cleaned = 0;
  for (size_t t = 0; t + 1 < feature->tags.size(); t += 2) {
    if (feature->tags[t] >= layer->keys.size()) {
      cerr << "Error: out of bounds feature key (" << feature->tags[t]
           << " in " << layer->keys.size() << endl;
      exit(EXIT_FAILURE);
    }
    if (feature->tags[t + 1] >= layer->values.size()) {
      cerr << "Error: out of bounds feature value (" << feature->tags[t + 1]
           << " in " << layer->values.size() << endl;
      exit(EXIT_FAILURE);
    }
    const string& key = layer->keys[feature->tags[t]];
    mvt_value* val = &layer->values[feature->tags[t + 1]];
    // Clear tags with value of type string and key having substring 'name'.
    if (val->type == mvt_string && !val->string_value.empty() &&
        key.find("name") != string::npos) {
      if (verbose) {
        cout << "Clearing [" << key.c_str() << "]:"
             << val->string_value.c_str() << endl;
      }
      val->string_value.clear();
      num_names_cleaned++;
    }
  }
  return num_names_cleaned;
}

// Returns the first feature of the first layer in the tile.
const mvt_feature* get_region_feature(const zxy& t,
                                      const mvt_tile& region_tile) {
  if (region_tile.layers.empty()) {
    cerr << "WARNING: Region tile " << t.z << " " << t.x << " " << t.y
         << " does not have any layer. Skipping clean-up of boundary tile."
         << endl;
    return NULL;
  } else if (region_tile.layers.size() > 1) {
    cerr << "WARNING: Region tile " << t.z << " " << t.x << " " << t.y
         << " has more than one layer: "
         << static_cast<int>(region_tile.layers.size()) << endl;
    return &region_tile.layers[0].features[0];
  } else if (region_tile.layers[0].features.empty()) {
    cerr << "WARNING: Region tile " << t.z << " " << t.x << " " << t.y
         << " does not have any feature. Skipping clean-up of boundary tile."
         << endl;
    return NULL;
  }
  return &region_tile.layers[0].features[0];
}

static const mvt_feature* region_feature = NULL;
static bool is_point_and_within_region(const mvt_feature& feature) {
  return region_feature != NULL && feature.type == VT_POINT &&
      feature_intersects_region(*region_feature, feature);
}

void clear_features_intersecting_region(
    const zxy& t, const mvt_tile& region_tile, mvt_tile* main_tile,
    int* num_point_features_removed, int* num_features_modified) {
  *num_point_features_removed = 0;
  *num_features_modified = 0;
  region_feature = get_region_feature(t, region_tile);

  // Loop through all layers and features in main_tile.
  // Apply the following rules:
  // 1. Remove any point feature within the region.
  // 2. Clear names on any feature intersecting the region.
  int num_names_cleaned = 0;
  for (size_t l = 0; l < main_tile->layers.size(); l++) {
    mvt_layer* layer = &main_tile->layers[l];
    // Remove point feature.
    vector<mvt_feature>* features = &layer->features;
    vector<mvt_feature>::const_iterator cit = features->end();
    features->erase(remove_if(features->begin(), features->end(),
                    is_point_and_within_region), features->end());
    *num_point_features_removed += cit - features->end();

    for (size_t f = 0; f < layer->features.size(); f++) {
      mvt_feature* feature = &layer->features[f];
      if (region_feature != NULL &&
          feature_intersects_region(*region_feature, *feature)) {
        if (feature->type != VT_POINT) {
          num_names_cleaned += clear_names(layer, feature);
          *num_features_modified += 1;
        }
      }
    }
  }
  if (verbose) {
    cout << "Modified tile " << t.z << " " << t.x << " " << t.y
         << ": removed " << *num_point_features_removed << " point features "
            "and cleaned " << num_names_cleaned << " names." << endl;
  }
}

void mbtiles_update_tile(sqlite3 *outdb, int tz, int tx, int ty,
    const char *data, int size) {
  sqlite3_stmt *stmt;
  const string query = "update tiles set tile_data = ? where zoom_level = ? "
      "and tile_column = ? and tile_row = ?;";
  if (sqlite3_prepare_v2(outdb, query.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
    cerr << "sqlite3 update prep failed." << endl;
    exit(EXIT_FAILURE);
  }
  sqlite3_bind_blob(stmt, 1, data, size, NULL);
  sqlite3_bind_int(stmt, 2, tz);
  sqlite3_bind_int(stmt, 3, tx);
  sqlite3_bind_int(stmt, 4, ty);
  if (sqlite3_step(stmt) != SQLITE_DONE) {
    cerr << "sqlite3 update failed: " << sqlite3_errmsg(outdb) << endl;
  }
  if (sqlite3_finalize(stmt) != SQLITE_OK) {
    cerr << "sqlite3 finalize failed: " << sqlite3_errmsg(outdb) << endl;
  }
}


void clear_boundary_tiles(const string& main_fname,
                          const string& region_fname,
                          set<zxy> boundary_tiles) {
  sqlite3 *main_db;
  sqlite3 *region_db;
  if (sqlite3_open(main_fname.c_str(), &main_db) != SQLITE_OK) {
    cerr << main_fname << ": " << sqlite3_errmsg(main_db) << endl;
    exit(EXIT_FAILURE);
  }
  if (sqlite3_open(region_fname.c_str(), &region_db) != SQLITE_OK) {
    cerr << region_fname << ": " << sqlite3_errmsg(region_db) << endl;
    exit(EXIT_FAILURE);
  }

  // Check all main features in boundary tiles and modify the ones
  // intersecting the region polygon.
  string compressed;
  string encoded;
  mvt_tile main_tile, region_tile;
  long long total_num_features_deleted = 0, total_num_features_modified = 0;
  const int kNumBoundaryTiles = boundary_tiles.size();
  int num_tiles_processed = 0;
  for (auto t : boundary_tiles) {
    read_tile(main_db, t, &main_tile);
    read_tile(region_db, t, &region_tile);
    int num_features_deleted = 0, num_features_modified = 0;
    clear_features_intersecting_region(t, region_tile, &main_tile,
        &num_features_deleted, &num_features_modified);
    total_num_features_deleted += num_features_deleted;
    total_num_features_modified += num_features_modified;

    // Write updated tile to db.
    encoded = main_tile.encode();
    compress(encoded, compressed);
    mbtiles_update_tile(main_db, t.z, t.x, t.y,
                        compressed.data(), compressed.size());
    num_tiles_processed++;
    if (num_tiles_processed == kNumBoundaryTiles) {
      cout << "All boundary tiles processed: " << num_tiles_processed << "/"
           << kNumBoundaryTiles << endl;
    } else if(num_tiles_processed%1000 == 0) {
      cout << "Boundary tiles processed: " << num_tiles_processed << "/"
           << kNumBoundaryTiles << "\r";
    }
  }
  // Close databases.
  if (sqlite3_close(main_db) != SQLITE_OK) {
    cerr << main_fname << ": could not close database: "
         << sqlite3_errmsg(main_db) << endl;
    exit(EXIT_FAILURE);
  }
  if (sqlite3_close(region_db) != SQLITE_OK) {
    cerr << region_fname << ": could not close database: "
         << sqlite3_errmsg(region_db) << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Total number of boundary point features deleted: "
       << total_num_features_deleted << endl;
  cout << "Total number of boundary features modified: "
       << total_num_features_modified << endl;
}

void mbtile_vacuum(const string& mbtile_name) {
  cerr << "Vacuuming sqlite3 db " << mbtile_name << "...";
  sqlite3* db;
  if (sqlite3_open(mbtile_name.c_str(), &db) != SQLITE_OK) {
     cerr << mbtile_name << ": " << sqlite3_errmsg(db) << endl;
     exit(EXIT_FAILURE);
  }
  char *err = 0;
  const string query = "vacuum;";
  int rc = sqlite3_exec(db, query.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << err << endl;
    sqlite3_free(err);
  }
  if (sqlite3_close(db) != SQLITE_OK) {
     cerr << mbtile_name << ": could not close database: "
          << sqlite3_errmsg(db) << endl;
     exit(EXIT_FAILURE);
  }
  cerr << " done." << endl;
}

void usage(char **argv) {
  cerr << "\nCuts out a region from an existing MBTiles (-m).\n"
          "The region to remove from the MBTiles is specified by two "
          "additional MBtiles (options -b and -r).\n\n"
          "Usage: "
       << argv[0]
       << " -m main.mbtiles -b boundary.mbtiles -r region.mbtiles "
          "[-h -R -B -V -v]\n\n"
          "Options:\n"
          "  --main-tileset\t\t-m\tMain MBtile to be modified (required).\n"
          "  --region-tileset\t\t-r\tRegion MBtile containing the region to "
          "cut out as polygon (required).\n"
          "  --boundary-tileset\t\t-b\tBoundary MBtile containing the outline "
          "of the region to cut out (required).\n"
          "  --help\t\t\t-h\tShow this usage.\n"
          "  --delete-tiles-within-region\t-R\tRemove tiles within region "
          "[default: off].\n"
          "  --change-boundary-features\t-B\tModify boundary tiles (clear "
          "names, remove interior features) [default: off].\n"
          "  --vacuum-mbtiles\t\t-V\tCompact modified MBtiles [default: off]."
       << endl << endl << endl;
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  extern char *optarg;
  int i;
  bool change_boundary_features = false;
  bool delete_tiles_within_region = false;
  bool vacuum_mbtiles = false;

  string main_tileset;
  string region_tileset;
  string boundary_tileset;
	
  struct option long_options[] = {
    {"main-tileset", required_argument, 0, 'm'},
	{"region-tileset", required_argument, 0, 'r'},
	{"boundary-tileset", required_argument, 0, 'b'},
	{"help", no_argument, 0, 'h'},
	{"change-boundary-features", no_argument, 0, 'B'},
	{"delete-tiles-within-region", no_argument, 0, 'R'},
	{"vacuum-mbtiles", no_argument, 0, 'V'},
	{"verbose", no_argument, 0, 'v'},
	{0, 0, 0, 0},
  };

  string getopt_str;
  for (size_t lo = 0; long_options[lo].name != NULL; lo++) {
	if (long_options[lo].val > ' ') {
      getopt_str.push_back(long_options[lo].val);
	  if (long_options[lo].has_arg == required_argument) {
        getopt_str.push_back(':');
	  }
	}
  }

  while ((i = getopt_long(argc, argv, getopt_str.c_str(), long_options,
         NULL)) != -1) {
    switch (i) {
	  case 0:
	    break;
	  case 'm':
	    main_tileset = optarg;
	    break;
      case 'r':
        region_tileset = optarg;
	    break;
      case 'b':
        boundary_tileset = optarg;
	    break;
      case 'h':
        usage(argv);
	    break;
      case 'v':
	    verbose = true;
	    break;
      case 'B':
        change_boundary_features = true;
	    break;
      case 'R':
        delete_tiles_within_region = true;
	    break;
      case 'V':
        vacuum_mbtiles = true;
	    break;
      default:
	    usage(argv);
	  }
	}

    if (main_tileset.empty() || region_tileset.empty() ||
        boundary_tileset.empty()) {
      usage(argv);
    }
    const chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    set<zxy> main_tiles;
    set<zxy> interior_tiles;
	set<zxy> boundary_tiles;	
    classify_tiles(
        region_tileset, boundary_tileset, main_tileset,
        &interior_tiles, &boundary_tiles, &main_tiles);

    if (delete_tiles_within_region) {
      delete_interior_tiles(main_tileset, main_tiles, interior_tiles);
    }
    if (change_boundary_features) {
      clear_boundary_tiles(main_tileset, region_tileset, boundary_tiles);
    }
    if (vacuum_mbtiles) {
      mbtile_vacuum(main_tileset);
    }

    const chrono::steady_clock::time_point end = chrono::steady_clock::now();
    const chrono::seconds secs =
        chrono::duration_cast<chrono::seconds>(end - begin);
    const unsigned int hours =
        static_cast<unsigned int>(secs.count() / 3600.0);
    cout << "Elapsed time: " << hours
         << " h " << static_cast<int>((secs.count() - 3600 * hours) / 60.0)
         << " m " << static_cast<int>(secs.count()%60) << " s" << endl;
	return 0;
}
