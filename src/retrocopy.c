#include "config.h"

#include <assert.h>
#include "db.h"
#include "log.h"
#include "retrocopy.h"

static sqlite3_stmt *
prepare_cluster_query_stmt (sqlite3 *db, const int filter)
{
	log_trace ("Inside %s", __func__);

	sqlite3_stmt *stmt = NULL;

	const char sql[] =
		"WITH\n"
		"	gene (gene_name, chr, start, end) AS (\n"
		"		SELECT gene_name, chr, MIN(start), MAX(end)\n"
		"		FROM exon\n"
		"		GROUP BY gene_name\n"
		"	),\n"
		"	gene_rank (gene_name, chr, start, end, dist) AS (\n"
		"		SELECT *,\n"
		"			DENSE_RANK() OVER (\n"
		"				PARTITION BY chr\n"
		"				ORDER BY start ASC, end ASC\n"
		"			)\n"
		"		FROM gene\n"
		"	)\n"
		"SELECT c.id, c.sid, c.chr, c.start, c.end,\n"
		"	c.gene_name, g.chr, g.start, g.end, g.dist\n"
		"FROM cluster AS c\n"
		"INNER JOIN gene_rank AS g\n"
		"	USING (gene_name)\n"
		"WHERE c.filter = $FILTER\n"
		"ORDER BY c.chr ASC, c.start ASC, c.end ASC";

	log_debug ("Query schema:\n%s", sql);
	stmt = db_prepare (db, sql);

	db_bind_int (stmt,
			sqlite3_bind_parameter_index (stmt, "$FILTER"),
			filter);

	return stmt;
}
