#include <getopt.h>
#include <sqlite3.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "geometry.hpp"
#include "mapbox/geometry/point.hpp"
#include "mapbox/geometry/polygon.hpp"
#include "mapbox/geometry/wagyu/util.hpp"
#include "mbtiles.hpp"
#include "mvt.hpp"
#include "protozero/pbf_reader.hpp"

using namespace std;
using mapbox::geometry::linear_ring;
using mapbox::geometry::point;
using mapbox::geometry::polygon;

bool verbose = false;

struct zxyrange {
  int zoom;
  long long col;
  long long rowmin;
  long long rowmax;
};

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
    //if (verbose) {
    //  cout << fname << " " << tile.z << " " << tile.x << " " <<  tile.y
    //       << " (flipped y: " << (1LL << tile.z) - 1 - tile.y << ")" << endl;
    //}
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
    set<zxy>* boundary_set, set<zxy>* outside_set, set<zxy>* main_set) {
  cerr << "Classifying tiles..." << endl;
  set<zxy> region_set;
  create_tile_set(region_fname, &region_set);
  create_tile_set(boundary_fname, boundary_set);
  set_difference(region_set.begin(), region_set.end(),
                 boundary_set->begin(), boundary_set->end(),
                 inserter(*interior_set, interior_set->end()));

  create_tile_set(main_fname, main_set);
  set_difference(main_set->begin(), main_set->end(),
                 region_set.begin(), region_set.end(),
                 inserter(*outside_set, outside_set->end()));

  const long long main = main_set->size();
  const long long interior = interior_set->size();
  const long long boundary = boundary_set->size();
  const long long outside = outside_set->size();
  cout << "Total number of tiles in main tileset: " << main << endl;
  cout << "Number of interior tiles: " << interior << endl;
  cout << "Number of boundary tiles: " << boundary << endl;
  cout << "Number of outside tiles: " << outside << endl;
  cout << "Difference (main-interior-boundary-outside): "
       << main - interior - boundary - outside << endl;
}

void identify_zoom_range(sqlite3* db, int* min_zoom, int* max_zoom) {
  if (db == NULL || min_zoom == NULL || max_zoom == NULL) {
    return;
  }
  string query = "select min(zoom_level) from tiles;";
  sqlite3_stmt* stmt;
  if (sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
    cerr << "select failed: " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  if (sqlite3_step(stmt) == SQLITE_ROW) {
    *min_zoom = sqlite3_column_int(stmt, 0);
  }
  sqlite3_finalize(stmt);

  query = "select max(zoom_level) from tiles;";
  if (sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
    cerr << "select failed: " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  if (sqlite3_step(stmt) == SQLITE_ROW) {
    *max_zoom = sqlite3_column_int(stmt, 0);
  }
  sqlite3_finalize(stmt);
  cout << "Zoom range of main tileset: " << *min_zoom << " - " << *max_zoom
       << endl;
}

void delete_ranges(zxyrange ranges, sqlite3* db, long long& num_tiles_deleted) {
  char query[200];
  snprintf(
      query, 200,
      "delete from tiles where zoom_level=%d and tile_column=%lld "
      "and (tile_row between %lld and %lld);",
      ranges.zoom, ranges.col, ranges.rowmin, ranges.rowmax);
  if (verbose) {
    cout << "[" << ranges.zoom << "/" << ranges.col << "/" << ranges.rowmin << "-" << ranges.rowmax << "]: "<< query
         << endl;
  }
  char *err = 0;
  int rc = sqlite3_exec(db, query, NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << err << endl;
    sqlite3_free(err);
  } else {
    num_tiles_deleted += (ranges.rowmax - ranges.rowmin + 1);
    cout << "Tiles deleted: " << num_tiles_deleted << '\r';
  }
   
};

void delete_tiles(string fname, set<zxy>* tiles_to_delete,
                           bool main_is_raster) {
  cerr << "Deleting tiles..." << endl;
  sqlite3 *db;
  if (sqlite3_open(fname.c_str(), &db) != SQLITE_OK) {
    cerr << fname << ": " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  // We assume that the size of raster tiles is 256 x 256 px. This results in
  // a zoom level increased by one and 2 x 2 tiles.
  const int num_tiles = main_is_raster ? 4 : 1;
  long long num_tiles_deleted = 0;
  int min_zoom = 0, max_zoom = 24;
  zxyrange ranges;
  ranges.zoom = -1;

  identify_zoom_range(db, &min_zoom, &max_zoom);
  for (auto t : *tiles_to_delete) {
    const int zoom_level = main_is_raster ? t.z + 1 : t.z;
    if (zoom_level < min_zoom || zoom_level > max_zoom) {
      continue;
    }
    for (int i = 0; i < num_tiles; ++i) {
      const long long tile_column = main_is_raster ? t.x * 2 + i%2 : t.x;
      const long long tile_row = main_is_raster ? t.y * 2 + i/2 : t.y;

      if (ranges.zoom == zoom_level &&
          ranges.col == tile_column &&
          ranges.rowmax == tile_row - 1) {
          ranges.rowmax = tile_row;
      } else {
        if (ranges.zoom != -1) {
          delete_ranges(ranges, db, num_tiles_deleted);
        }
        ranges.zoom = zoom_level;
        ranges.col = tile_column;
        ranges.rowmin = tile_row;
        ranges.rowmax = tile_row;
      }
    }
  }
  if (ranges.zoom != -1) {
    delete_ranges(ranges, db, num_tiles_deleted);
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

void get_region_polygons(const zxy& t,
                         const mvt_tile& region_tile,
                         vector<polygon<int64_t> >* polygons) {
  if (polygons == NULL) {
    return;
  }
  polygons->clear();
  if (region_tile.layers.empty()) {
    cerr << "WARNING: Region tile " << t.z << " " << t.x << " " << t.y
         << " does not have any layer. Skipping clean-up of boundary tile."
         << endl;
    return;
  } else if (region_tile.layers[0].features.empty()) {
    cerr << "WARNING: Region tile " << t.z << " " << t.x << " " << t.y
         << " does not have any feature. Skipping clean-up of boundary tile."
         << endl;
    return;
  } else if (region_tile.layers.size() > 1) {
    cerr << "WARNING: Region tile " << t.z << " " << t.x << " " << t.y
         << " has more than one layer: "
         << static_cast<int>(region_tile.layers.size()) << endl;
  }
  if (region_tile.layers[0].features[0].type != mvt_polygon) {
    cerr << "WARNING: Feature in region tile " << t.z << " " << t.x << " "
         << t.y << " does not have polygon geometry. Skipping clean-up of "
         "boundary tile." << endl;
  }
  int64_t num_rings = 0;
  for (size_t f = 0; f < region_tile.layers[0].features.size(); f++) {
    polygon<int64_t> poly;
    linear_ring<int64_t> ring;
    const mvt_feature& feature = region_tile.layers[0].features[f];
    if (feature.type != mvt_polygon) {
      continue;
    }
    for (size_t g = 0; g < feature.geometry.size(); g++) {
      if (feature.geometry[g].op == VT_MOVETO) {
        // Start of new ring.
        ring.push_back({feature.geometry[g].x, feature.geometry[g].y});
      } else if (feature.geometry[g].op == VT_CLOSEPATH) {
        // Close ring by adding first point.
        ring.push_back(ring.front());
        // Add recent ring.
        poly.push_back(ring);
        ring.clear();
        num_rings++;
      } else {
        ring.push_back({feature.geometry[g].x, feature.geometry[g].y});
      }
    }
    polygons->push_back(poly);
  }
  if (verbose) {
    cerr << "Tile [" << t.z << " " << t.x << " " << t.y
         << "]: num_polygons: " << polygons->size() << ", num_rings: "
         << num_rings << endl;
  }
}

// Modified from https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
static int pnpoly(const linear_ring<int64_t>& ring,
                  const point<int64_t>& point) {
  size_t i, j;
  bool c = false;
  for (i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
    if (((ring[i].y > point.y) != (ring[j].y > point.y)) &&
        (point.x < (ring[j].x - ring[i].x) * (point.y - ring[i].y) /
         (double) (ring[j].y - ring[i].y) + ring[i].x))  {
      c = !c;
    }
  }
  return c;
}

bool feature_inside_region(vector<polygon<int64_t> > region_polygons,
                               const mvt_feature& main_feature) {
  if (region_polygons.empty()) {
    return false;
  }
  // A feature is inside if all points lie within a polygon.
  // This is an approximation which should be good enough for our purpose.
  for (size_t g = 0; g < main_feature.geometry.size(); g++) {
    const point<int64_t> pt(main_feature.geometry[g].x,
                            main_feature.geometry[g].y);
    // Check the point against every ring in polygons.
    for (size_t p = 0; p < region_polygons.size(); p++) {
      for (size_t r = 0; r < region_polygons[p].size(); r++) {
        const linear_ring<int64_t>& ring = region_polygons[p][r];
        // TODO: Add support for holes where a ring could be a inner ring,
        // i.e. a hole.
        // const bool is_outer_ring = area(ring) > 0 ? true : false;
        //
        // return false as soon as a point outside is found.
        if (!pnpoly(ring, pt)) {
          return false;
        }
      }
    }
  }
  return true;
}

bool feature_intersects_region(vector<polygon<int64_t> > region_polygons,
                               const mvt_feature& main_feature) {
  if (region_polygons.empty()) {
    return false;
  }
  // A feature intersects as soon as one point lies within a polygon.
  for (size_t g = 0; g < main_feature.geometry.size(); g++) {
    const point<int64_t> pt(main_feature.geometry[g].x,
                            main_feature.geometry[g].y);
    // Check the point against every ring in polygons.
    for (size_t p = 0; p < region_polygons.size(); p++) {
      for (size_t r = 0; r < region_polygons[p].size(); r++) {
        const linear_ring<int64_t>& ring = region_polygons[p][r];
        // TODO: Add support for holes where a ring could be a inner ring,
        // i.e. a hole.
        // const bool is_outer_ring = area(ring) > 0 ? true : false;
        if (pnpoly(ring, pt)) {
          return true;
        }
      }
    }
  }
  return false;
}

static vector<polygon<int64_t>> region_polygons;
static bool is_within_region(const mvt_feature& feature) {
  return !region_polygons.empty() &&
      feature_inside_region(region_polygons, feature);
}

void clear_features_intersecting_region(
    const zxy& t, const mvt_tile& region_tile, mvt_tile* main_tile,
    int* num_point_features_removed, int* num_features_modified) {
  *num_point_features_removed = 0;
  *num_features_modified = 0;
  get_region_polygons(t, region_tile, &region_polygons);

  // Loop through all layers and features in main_tile.
  // Apply the following rules:
  // 1. Remove any point feature within the region.
  // 2. Remove any feature fully within the region.
  // 3. Clear names on any feature intersecting the region.
  int num_names_cleaned = 0, num_point_features_kept = 0,
    num_lines_polygons_not_modified = 0;
  for (size_t l = 0; l < main_tile->layers.size(); l++) {
    mvt_layer* layer = &main_tile->layers[l];
    // Remove point feature.
    vector<mvt_feature>* features = &layer->features;
    vector<mvt_feature>::const_iterator cit = features->end();
    features->erase(remove_if(features->begin(), features->end(),
                    is_within_region), features->end());
    *num_point_features_removed += cit - features->end();

    // Clear names on polyline or polygon features.
    for (size_t f = 0; f < layer->features.size(); f++) {
      mvt_feature* feature = &layer->features[f];
      if (feature->type == VT_POINT) {
        num_point_features_kept++;
      } else if (feature_intersects_region(region_polygons, *feature)) {
        num_names_cleaned += clear_names(layer, feature);
        *num_features_modified += 1;
      } else {
        num_lines_polygons_not_modified++;
      }
    }
  }
  if (verbose) {
    cout << "Modified tile " << t.z << " " << t.x << " " << t.y
         << ": Point features removed: " << *num_point_features_removed
         << ", kept: " << num_point_features_kept
         << " Line/Polygon features modified: " << *num_features_modified
         << ", not modified: " << num_lines_polygons_not_modified
         << ", names cleaned: " << num_names_cleaned << endl;
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
  cerr << "Process features in boundary tiles..." << endl;
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
  cerr << "Vacuum sqlite3 db " << mbtile_name << "...";
  sqlite3* db;
  if (sqlite3_open(mbtile_name.c_str(), &db) != SQLITE_OK) {
     cerr << mbtile_name << ": " << sqlite3_errmsg(db) << endl;
     exit(EXIT_FAILURE);
  }
  char* err = 0;
  const string kQuery("vacuum;");
  int rc = sqlite3_exec(db, kQuery.c_str(), NULL, NULL, &err);
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

void clone_mbtile(const string& original_mbtile, const string& main_mbtile) {
  cerr << "Copying tiles & metadata tables from " << original_mbtile
       << " into " << main_mbtile << "...\n";
  sqlite3* db;
  if (sqlite3_open(main_mbtile.c_str(), &db) != SQLITE_OK) {
    cerr << main_mbtile << ": " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  char* err = 0;
  // Attach original MBTile as originaldb.
  const string kAttachQuery("ATTACH DATABASE \"" + original_mbtile +
      "\" AS originaldb;");
  int rc = sqlite3_exec(db, kAttachQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kAttachQuery << endl << err << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }
  // Create tiles table.
  const string kCreateTilesQuery("CREATE TABLE tiles (zoom_level integer,"
      "tile_column integer, tile_row integer, tile_data blob);");
  rc = sqlite3_exec(db, kCreateTilesQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kCreateTilesQuery << endl << err << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }
  // Copy tiles table from originaldb.
  const string kCopyTilesQuery("INSERT INTO main.tiles(zoom_level, "
      "tile_column, tile_row, tile_data) SELECT zoom_level, tile_column, "
      "tile_row, tile_data FROM originaldb.tiles;");
  rc = sqlite3_exec(db, kCopyTilesQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kCopyTilesQuery << endl << err << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }
  // Create tiles index.
  const string kIndexTilesQuery("CREATE UNIQUE INDEX tile_index on "
      "tiles (zoom_level, tile_column, tile_row);");
  rc = sqlite3_exec(db, kIndexTilesQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kIndexTilesQuery << endl << err << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }

  // Create metadata table.
  const string kCreateMetadataQuery(
      "CREATE TABLE metadata (name text, value text);");
  rc = sqlite3_exec(db, kCreateMetadataQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kCreateMetadataQuery << endl << err << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }
  // Copy metadata table from originaldb.
  const string kCopyMetadataQuery("INSERT INTO main.metadata (name, value) "
      "SELECT name, value FROM originaldb.metadata;");
  rc = sqlite3_exec(db, kCopyMetadataQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kCopyMetadataQuery << endl << err << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }
  // Create metadata index.
  const string kIndexMetadataQuery(
      "CREATE UNIQUE INDEX name on metadata (name);");
  rc = sqlite3_exec(db, kIndexMetadataQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << kIndexMetadataQuery << endl << cerr << endl;
    sqlite3_free(err);
    exit(EXIT_FAILURE);
  }

  // Close main MBTile.
  if (sqlite3_close(db) != SQLITE_OK) {
     cerr << main_mbtile << ": could not close database: "
          << sqlite3_errmsg(db) << endl;
     exit(EXIT_FAILURE);
  }
}

void update_tiles_index(const string& main_mbtile) {
  cerr << "Updating tiles index..." << endl;
  sqlite3* db;
  if (sqlite3_open(main_mbtile.c_str(), &db) != SQLITE_OK) {
    cerr << main_mbtile << ": " << sqlite3_errmsg(db) << endl;
    exit(EXIT_FAILURE);
  }
  char* err = 0;
  // Create tiles index.
  const string kIndexTilesQuery("REINDEX tile_index;");
  int rc = sqlite3_exec(db, kIndexTilesQuery.c_str(), NULL, NULL, &err);
  if (rc != SQLITE_OK){
    cerr << "SQL error: " << err << endl;
    sqlite3_free(err);
  }

  // Close main MBTile.
  if (sqlite3_close(db) != SQLITE_OK) {
     cerr << main_mbtile << ": could not close database: "
          << sqlite3_errmsg(db) << endl;
     exit(EXIT_FAILURE);
  }
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
          "  --original-tileset\t\t-o\tOptional MBtile as input. If set, this "
          "tileset acts as source for the newly created main MBtile.\n"
          "  --main-is-raster\t\t-s\tIf set the main tileset contains raster "
          "data with a tile size of 256 px. [default: off]\n"
          "  --help\t\t\t-h\tShow this usage.\n"
          "  --delete-tiles-within-region\t-R\tRemove tiles within region "
          "[default: off].\n"
          "  --change-boundary-features\t-B\tModify boundary tiles (clear "
          "names, remove interior features) [default: off].\n"
          "  --vacuum-mbtiles\t\t-V\tCompact modified MBtiles [default: off]."
       << endl << endl << endl;
  exit(EXIT_FAILURE);
}

void print_elapsed_time(const chrono::steady_clock::time_point& begin) {
  const chrono::steady_clock::time_point end = chrono::steady_clock::now();
  const chrono::seconds secs =
      chrono::duration_cast<chrono::seconds>(end - begin);
  const unsigned int hours =
      static_cast<unsigned int>(secs.count() / 3600.0);
  cerr << "Elapsed time: " << hours
       << " h " << static_cast<int>((secs.count() - 3600 * hours) / 60.0)
       << " m " << static_cast<int>(secs.count()%60) << " s" << endl << endl;
}

int main(int argc, char **argv) {
  extern char *optarg;
  int i;
  bool change_boundary_features = false;
  bool delete_tiles_within_region = false;
  bool delete_tiles_outside_region = false;
  bool main_is_raster = false;
  bool vacuum_mbtiles = false;

  string main_tileset_fname;
  string region_tileset_fname;
  string boundary_tileset_fname;
  string original_tileset_fname;
	
  struct option long_options[] = {
    {"main-tileset", required_argument, 0, 'm'},
	{"region-tileset", required_argument, 0, 'r'},
	{"boundary-tileset", required_argument, 0, 'b'},
	{"original-tileset", required_argument, 0, 'o'},
	{"main-is-raster", no_argument, 0, 's'},
	{"help", no_argument, 0, 'h'},
	{"change-boundary-features", no_argument, 0, 'B'},
	{"delete-tiles-within-region", no_argument, 0, 'R'},
	{"delete-tiles-outside-boundary", no_argument, 0, 'O'},
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
        main_tileset_fname = optarg;
        break;
      case 'r':
        region_tileset_fname = optarg;
        break;
      case 'b':
        boundary_tileset_fname = optarg;
        break;
      case 'o':
        original_tileset_fname = optarg;
        break;
      case 's':
        main_is_raster = true;
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
      case 'O':
        delete_tiles_outside_region = true;
      case 'V':
        vacuum_mbtiles = true;
        break;
      default:
        usage(argv);
      }
  }

  if (main_tileset_fname.empty() || region_tileset_fname.empty() ||
      boundary_tileset_fname.empty()) {
    usage(argv);
  }
  const chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  if (!original_tileset_fname.empty()) {
    clone_mbtile(original_tileset_fname, main_tileset_fname);
    print_elapsed_time(begin);
  }

  set<zxy> main_tiles;
  set<zxy> interior_tiles;
  set<zxy> boundary_tiles;
  set<zxy> outside_tiles;
  classify_tiles(
      region_tileset_fname, boundary_tileset_fname, main_tileset_fname,
      &interior_tiles, &boundary_tiles, &outside_tiles, &main_tiles);
  print_elapsed_time(begin);
  cerr << "Classifying done..." << endl;

  if (delete_tiles_outside_region) {
    cerr << "Deleting tiles outside region..." << endl;
    delete_tiles(main_tileset_fname, &outside_tiles, main_is_raster);
    print_elapsed_time(begin);
    cerr << "Deleted tiles outside region..." << endl;
  }

  if (delete_tiles_within_region) {
    cerr << "Deleting tiles inside region..." << endl;
    delete_tiles(main_tileset_fname, &interior_tiles, main_is_raster);
    print_elapsed_time(begin);
    cerr << "Deletied tiles inside region..." << endl;
  }
  cerr << "Deletions over (if any)..." << endl;

  if (change_boundary_features) {
    cerr << "Changing boundary features..." << endl;
    if (main_is_raster) {
      cerr << "WARNING: Processing of raster boundary tiles skipped." << endl;
    } else {
      cerr << "clear boundary tiles..." << endl;
      clear_boundary_tiles(main_tileset_fname, region_tileset_fname,
                           boundary_tiles);
      print_elapsed_time(begin);
      cerr << "cleared boundary tiles..." << endl;
    }
  }
  // Finalize MBtile.
  cerr << "Finalising..." << endl;
  if (!original_tileset_fname.empty() &&
      (delete_tiles_within_region || delete_tiles_outside_region || change_boundary_features)) {
    update_tiles_index(main_tileset_fname);
    print_elapsed_time(begin);
  }
  if (vacuum_mbtiles) {
    mbtile_vacuum(main_tileset_fname);
    print_elapsed_time(begin);
  }
  return 0;
}
