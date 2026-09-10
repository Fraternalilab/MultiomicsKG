_examples = [
    {
        "question": "List the pairs of proteins that are part of a protein complex.",
        "query": """
            MATCH (p1:protein)-[r:interacts]->(p2:protein)
            MATCH (p1)-[:maps]->(pdb:pdb)<-[:maps]-(p2)
            RETURN p1.name AS Protein1, p2.name AS Protein2, pdb.id AS PDB_ID
            LIMIT 10
        """,
    },
    {
        "question": "A protein complex is a pair of proteins that interact and both map to the same pdb. List the protein complexes.",
        "query": """
            MATCH (p1:protein)-[r:interacts]->(p2:protein)
            MATCH (p1)-[:maps]->(pdb:pdb)<-[:maps]-(p2)
            RETURN p1.name AS Protein1, p2.name AS Protein2, pdb.id AS PDB_ID
            LIMIT 10
        """,
    },
]


def get_fewshot_examples():
    """
    Returns a list of few-shot examples for the LangGraph workflow.
    """
    return _examples
