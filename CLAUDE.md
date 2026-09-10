## Neo4j database schema

Before generating Cypher queries:

1. Read neo4j-schema.json.
2. Use the cached schema to identify node labels,
   relationship types, and property names.
3. Do not call the Neo4j get-schema tool unless the
   cached schema is missing or needs refreshing.
4. Refresh the cached schema after database schema changes.