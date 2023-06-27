/*
 * sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
 * Copyright (C) 2019-2023 Thiago L. A. Miller <tmiller@mochsl.org.br
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef KV_H
#define KV_H

typedef void (*KVFunc) (const char *key, const void *value, void *user_data);

typedef struct _KV KV;

KV         * kv_new       (const char *path);
void         kv_free      (KV *kv);

const char * kv_path      (KV *kv);
int          kv_count     (KV *kv);

void         kv_insert    (KV *kv, const char *key, const void *value, int n);
void         kv_del_key   (KV *kv, const char *key);
const void * kv_get_value (KV *kv, const char *key);

void         kv_foreach   (KV *kv, KVFunc func, void *user_data);

#define kv_contains(kv,key) (kv_get_value ((kv),(key)) != NULL) 

#endif /* kv.h */
