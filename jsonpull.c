#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include "jsonpull.h"

typedef enum json_expect {
	JSON_ITEM, JSON_COMMA, JSON_COLON, JSON_KEY, JSON_VALUE,
} json_expect;

static int read_file(json_pull *p) {
	return fgetc(p->source);
}

static int peek_file(json_pull *p) {
	int c = getc(p->source);
	ungetc(c, p->source);
	return c;
}

json_pull *json_begin_file(FILE *f) {
	json_pull *j = malloc(sizeof(json_pull));
	j->container = NULL;

	j->read = read_file;
	j->peek = peek_file;
	j->source = f;

	return j;
}

static int read_string(json_pull *p) {
	char *cp = p->source;
	if (*cp == '\0') {
		return EOF;
	}
	int c = (unsigned char) *cp;
	cp++;
	p->source = cp;
	return c;
}

static int peek_string(json_pull *p) {
	char *cp = p->source;
	if (*cp == '\0') {
		return EOF;
	}
	return (unsigned char) *cp;
}

json_pull *json_begin_string(char *s) {
	json_pull *j = malloc(sizeof(json_pull));
	j->container = NULL;

	j->read = read_string;
	j->peek = peek_string;
	j->source = s;

	return j;
}

#define SIZE_FOR(i) (((i) + 31) & ~31)

static json_object *add_object(json_pull *j, json_type type) {
	json_object *o = malloc(sizeof(struct json_object));
	o->type = type;
	o->parent = j->container;
	o->array = NULL;
	o->keys = NULL;
	o->values = NULL;
	o->length = 0;

	json_object *c = j->container;

	if (c != NULL) {
		if (c->type == JSON_ARRAY) {
			if (SIZE_FOR(c->length + 1) != SIZE_FOR(c->length)) {
				c->array = realloc(c->array, SIZE_FOR(c->length + 1) * sizeof(json_object *));
			}

			c->array[c->length++] = o;
			c->expect = JSON_COMMA;
		} else if (c->type == JSON_HASH) {
			if (c->expect == JSON_VALUE) {
				c->values[c->length - 1] = o;
				c->expect = JSON_COMMA;
			} else {
				if (type != JSON_STRING) {
					j->error = "Hash key is not a string";
				}

				if (SIZE_FOR(c->length + 1) != SIZE_FOR(c->length)) {
					c->keys = realloc(c->keys, SIZE_FOR(c->length + 1) * sizeof(json_object *));
					c->values = realloc(c->values, SIZE_FOR(c->length + 1) * sizeof(json_object *));
				}

				c->keys[c->length] = o;
				c->values[c->length] = NULL;
				c->length++;
				c->expect = JSON_COLON;
			}
		}
	}

	return o;
}

json_object *json_hash_get(json_object *o, char *s) {
	if (o == NULL || o->type != JSON_HASH) {
		return NULL;
	}

	int i;
	for (i = 0; i < o->length; i++) {
		if (o->keys[i] != NULL && o->keys[i]->type == JSON_STRING) {
			if (strcmp(o->keys[i]->string, s) == 0) {
				return o->values[i];
			}
		}
	}

	return NULL;
}

struct string {
	char *buf;
	int n;
	int nalloc;
};

static void string_init(struct string *s) {
	s->nalloc = 500;
	s->buf = malloc(s->nalloc);
	s->n = 0;
	s->buf[0] = '\0';
}

static void string_append(struct string *s, char c) {
	if (s->n + 2 >= s->nalloc) {
		s->nalloc += 500;
		s->buf = realloc(s->buf, s->nalloc);
	}

	s->buf[s->n++] = c;
	s->buf[s->n] = '\0';
}

static void string_free(struct string *s) {
	free(s->buf);
}

json_object *json_parse(json_pull *j) {
	int c;
again:
	/////////////////////////// Whitespace

	do {
		c = j->read(j);
		if (c == EOF) {
			if (j->container != NULL) {
				j->error = "Reached EOF without all containers being closed";
			}

			return NULL;
		}
	} while (c == ' ' || c == '\t' || c == '\r' || c == '\n');

	/////////////////////////// Arrays

	if (c == '[') {
		j->container = add_object(j, JSON_ARRAY);
		j->container->expect = JSON_ITEM;
		goto again;
	} else if (c == ']') {
		if (j->container == NULL) {
			j->error = "Found ] at top level";
			return NULL;
		}

		if (j->container->type != JSON_ARRAY) {
			j->error = "Found ] not in an array";
			return NULL;
		}

		if (j->container->expect != JSON_COMMA) {
			if (! (j->container->expect == JSON_ITEM && j->container->length == 0)) {
				j->error = "Found ] without final element";
				return NULL;
			}
		}

		json_object *ret = j->container;
		j->container = ret->parent;
		return ret;
	}

	/////////////////////////// Hashes

	if (c == '{') {
		j->container = add_object(j, JSON_HASH);
		j->container->expect = JSON_KEY;
		goto again;
	} else if (c == '}') {
		if (j->container == NULL) {
			j->error = "Found } at top level";
			return NULL;
		}

		if (j->container->type != JSON_HASH) {
			j->error = "Found } not in a hash";
			return NULL;
		}

		if (j->container->expect != JSON_COMMA) {
			if (! (j->container->expect == JSON_KEY && j->container->length == 0)) {
				j->error = "Found } without final element";
				return NULL;
			}
		}

		json_object *ret = j->container;
		j->container = ret->parent;
		return ret;
	}

	/////////////////////////// Null

	if (c == 'n') {
		if (j->read(j) != 'u' || j->read(j) != 'l' || j->read(j) != 'l') {
			j->error = "Found misspelling of null";
			return NULL;
		}

		return add_object(j, JSON_NULL);
	}

	/////////////////////////// True

	if (c == 't') {
		if (j->read(j) != 'r' || j->read(j) != 'u' || j->read(j) != 'e') {
			j->error = "Found misspelling of true";
			return NULL;
		}

		return add_object(j, JSON_TRUE);
	}

	/////////////////////////// False

	if (c == 'f') {
		if (j->read(j) != 'a' || j->read(j) != 'l' || j->read(j) != 's' || j->read(j) != 'e') {
			j->error = "Found misspelling of false";
			return NULL;
		}

		return add_object(j, JSON_FALSE);
	}

	/////////////////////////// Comma

	if (c == ',') {
		if (j->container == NULL) {
			j->error = "Found comma at top level";
			return NULL;
		}

		if (j->container->expect != JSON_COMMA) {
			j->error = "Found unexpected comma";
			return NULL;
		}

		if (j->container->type == JSON_HASH) {
			j->container->expect = JSON_KEY;
		} else {
			j->container->expect = JSON_ITEM;
		}

		goto again;
	}

	/////////////////////////// Colon

	if (c == ':') {
		if (j->container == NULL) {
			j->error = "Found colon at top level";
			return NULL;
		}

		if (j->container->expect != JSON_COLON) {
			j->error = "Found unexpected colon";
			return NULL;
		}

		j->container->expect = JSON_VALUE;
		goto again;
	}

	/////////////////////////// Numbers

	if (c == '-' || (c >= '0' && c <= '9')) {
		struct string val;
		string_init(&val);

		if (c == '-') {
			string_append(&val, c);
			c = j->read(j);
		}

		if (c == '0') {
			string_append(&val, c);
		} else if (c >= '1' && c <= '9') {
			string_append(&val, c);
			c = j->peek(j);

			while (c >= '0' && c <= '9') {
				string_append(&val, j->read(j));
				c = j->peek(j);
			}
		}

		if (j->peek(j) == '.') {
			string_append(&val, j->read(j));

			c = j->peek(j);
			while (c >= '0' && c <= '9') {
				string_append(&val, j->read(j));
				c = j->peek(j);
			}
		}

		c = j->peek(j);
		if (c == 'e' || c == 'E') {
			string_append(&val, j->read(j));

			c = j->peek(j);
			if (c == '+' || c == '-') {
				string_append(&val, j->read(j));
			}

			c = j->peek(j);
			if (c < '0' || c > '9') {
				j->error = "Exponent without digits";
				return NULL;
			}
			while (c >= '0' && c <= '9') {
				string_append(&val, j->read(j));
				c = j->peek(j);
			}
		}

		json_object *n = add_object(j, JSON_NUMBER);
		n->number = atof(val.buf);
		string_free(&val);
		return n;
	}

	/////////////////////////// Strings

	if (c == '"') {
		struct string val;
		string_init(&val);

		while ((c = j->read(j)) != EOF) {
			if (c == '"') {
				break;
			} else if (c == '\\') {
				c = j->read(j);

				if (c == '"') {
					string_append(&val, '"');
				} else if (c == '\\') {
					string_append(&val, '\\');
				} else if (c == '/') {
					string_append(&val, '/');
				} else if (c == 'b') {
					string_append(&val, '\b');
				} else if (c == 'f') {
					string_append(&val, '\f');
				} else if (c == 'n') {
					string_append(&val, '\n');
				} else if (c == 'r') {
					string_append(&val, '\r');
				} else if (c == 't') {
					string_append(&val, '\t');
				} else if (c == 'u') {
					char hex[5] = "aaaa";
					int i;
					for (i = 0; i < 4; i++) {
						hex[i] = j->read(j);
					}
					unsigned long ch = strtoul(hex, NULL, 16);
					if (ch <= 0x7F) {
						string_append(&val, ch);
					} else if (ch <= 0x7FF) {
						string_append(&val, 0xC0 | (ch >> 6));
						string_append(&val, 0x80 | (ch & 0x3F));
					} else {
						string_append(&val, 0xE0 | (ch >> 12));
						string_append(&val, 0x80 | ((ch >> 6) & 0x3F));
						string_append(&val, 0x80 | (ch & 0x3F));
					}
				} else {
					j->error = "Found backslash followed by unknown character";
					return NULL;
				}
			} else {
				string_append(&val, c);
			}
		}

		json_object *s = add_object(j, JSON_STRING);
		s->string = val.buf;
		s->length = val.n;
		return s;
	}

	j->error = "Found unexpected character";
	return NULL;
}

void json_free(json_object *o) {
	int i;

	if (o == NULL) {
		return;
	}

	// Free any data linked from here

	if (o->type == JSON_ARRAY) {
		for (i = 0; i < o->length; i++) {
			json_free(o->array[i]);
		}
	} else if (o->type == JSON_HASH) {
		for (i = 0; i < o->length; i++) {
			json_free(o->keys[i]);
			json_free(o->values[i]);
		}
	} else if (o->type == JSON_STRING) {
		free(o->string);
	}

	// Expunge references to this as an array element
	// or a hash key or value.

	if (o->parent != NULL) {
		if (o->parent->type == JSON_ARRAY) {
			for (i = 0; i < o->parent->length; i++) {
				if (o->parent->array[i] == o) {
					break;
				}
			}

			if (i < o->parent->length) {
				memmove(o->parent->array + i, o->parent->array + i + 1, o->parent->length - i - 1);
				o->parent->length--;
			}
		}

		if (o->parent->type == JSON_HASH) {
			for (i = 0; i < o->parent->length; i++) {
				if (o->parent->keys[i] == o || o->parent->values[i] == o) {
					break;
				}
			}

			if (i < o->parent->length) {
				json_object *k = o->parent->keys[i];
				json_object *v = o->parent->values[i];

				memmove(o->parent->keys + i, o->parent->keys + i + 1, o->parent->length - i - 1);
				memmove(o->parent->values + i, o->parent->values + i + 1, o->parent->length - i - 1);
				o->parent->length--;
				
				if (o->parent->values[i] == NULL) {
					// Have not yet read the value for this pair
				} else {
					if (o->parent->keys[i] == o) {
					} else {
					}
				}
			}
		}
	}

	free(o);
}