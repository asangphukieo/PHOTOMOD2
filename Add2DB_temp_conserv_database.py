import sqlite3
conn = sqlite3.connect("GO_neighbor_enrichment.db")
c = conn.cursor()
#c.execute("DROP TABLE neighbor_conservation")
c.execute("CREATE TABLE neighbor_conservation (Query_cluster_id INTEGER NOT NULL, neighbor_cluster_id INTEGER NOT NULL, neighbor_found INTEGER NOT NULL, qeury_found INTEGER NOT NULL, number_of_genome_found INTEGER NOT NULL, C_score DOUBLE NOT NULL, D_score DOUBLE NOT NULL, CD_score DOUBLE NOT NULL , FOREIGN KEY(Query_cluster_id) REFERENCES GO_neighbor_enrichment (Query_cluster_id), FOREIGN KEY(neighbor_cluster_id) REFERENCES GO_neighbor_enrichment (neighbor_cluster_id), PRIMARY KEY (Query_cluster_id, neighbor_cluster_id) )")


db_list = []
with open("./CD_score_Qid_all.out") as infile:
    for line in infile:
        data = line.replace('\n','')
	data = data.split()
        db_list.append(data)

c.executemany("INSERT INTO neighbor_conservation VALUES (?, ?, ?, ?, ?, ?, ?, ?)", db_list)
c.execute("CREATE INDEX idx4 ON neighbor_conservation(Query_cluster_id, neighbor_cluster_id)")
conn.commit()
conn.close()
