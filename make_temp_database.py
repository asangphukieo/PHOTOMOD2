import sqlite3
conn = sqlite3.connect("GO_neighbor_enrichment.db")
c = conn.cursor()
c.execute("CREATE TABLE GO_neighbor_enrichment (Query_cluster_id INTEGER NOT NULL, neighbor_cluster_id INTEGER NOT NULL, total_neighbor_cluster_id INTEGER NOT NULL, x INTEGER NOT NULL, M INTEGER NOT NULL, K INTEGER NOT NULL, N INTEGER NOT NULL, pvalue DOUBLE NOT NULL, pvalue_adjust DOUBLE NOT NULL)")

db_list = []
with open("./neighbor_enrichment_no_header.out") as infile:
    for line in infile:
        data = line.replace('\n','')
	data = data.split()
        db_list.append(data)

c.executemany("INSERT INTO GO_neighbor_enrichment VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", db_list)
c.execute("CREATE INDEX idx1 ON GO_neighbor_enrichment(Query_cluster_id,neighbor_cluster_id)")
conn.commit()
conn.close()
