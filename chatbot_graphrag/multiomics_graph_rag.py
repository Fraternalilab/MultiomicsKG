import os
import logging

from langchain_neo4j import Neo4jGraph
from langchain_openai import ChatOpenAI
from langchain_ollama import ChatOllama
from langchain_anthropic import ChatAnthropic
from chatbot_graphrag.mgr.rag_workflow import setup_langgraph
from pydantic import SecretStr

RECURSION_LIMIT = (
    100  # Set a recursion limit for the LangGraph workflow to prevent infinite loops.
)
logging.basicConfig(level=logging.ERROR)


class MultiomicsGraphRAG:
    def __init__(
        self,
        enhanced_graph: Neo4jGraph,
        refresh_schema: bool,
        provider: str,
        model: str,
        api_key: str,
    ):
        if provider == "openai":
            self.llm = ChatOpenAI(
                model=model, temperature=0, api_key=SecretStr(api_key)
            )
        elif provider == "ollama":
            self.llm = ChatOllama(model=model, temperature=0)
        elif provider == "anthropic":
            self.llm = ChatAnthropic(model_name=model, timeout=120, stop=None, api_key=SecretStr(api_key))
        else:
            raise ValueError(
                f"Unsupported provider: {provider}. Supported providers are 'openai', 'ollama' and 'anthropic'."
            )
        self.langgraph = setup_langgraph(
            llm=self.llm,
            enhanced_graph=enhanced_graph,
            refresh_schema=refresh_schema,
            provider=provider,
        )

    def query(self, question: str) -> dict:
        """
        Process a question through the LangGraph workflow and return the final answer.
        """
        input_state = {"question": question}

        output_state = self.langgraph.invoke(
            input_state, {"recursion_limit": RECURSION_LIMIT}
        )
        return output_state


if __name__ == "__main__":
    NEO4J_URI = os.getenv("NEO4J_URI", "bolt://neo4j:7687")
    NEO4J_USERNAME = os.getenv("NEO4J_USERNAME", "neo4j")
    NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "snvpredsnvpred")
    NEO4J_DATABASE = os.getenv("NEO4J_DATABASE", "neo4j")
    OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o-mini")
    OLLAMA_MODEL = os.getenv("OLLAMA_MODEL", "gemma3:27b")
    OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "")

    enhanced_graph = Neo4jGraph(
        url=NEO4J_URI,
        username=NEO4J_USERNAME,
        password=NEO4J_PASSWORD,
        database=NEO4J_DATABASE,
        enhanced_schema=True,
        timeout=120,
    )

    rag = MultiomicsGraphRAG(
        enhanced_graph=enhanced_graph,
        refresh_schema=True,
        provider="openai",
        model=OPENAI_MODEL,
        api_key=OPENAI_API_KEY,
    )

    # question = "What's the weather in the UK?"
    # result = rag.query(question)
    # print("Query Result:", result)

    question = (
        "Find me proteins that interact with their interact relationship property of string combined score greater than 800, "
        "then find me the pair of proteins that map to the same pdb. "
        "Show me the protein names and the pdb ids in a table with up to 10 items."
    )
    result = rag.query(question)
    print("Query Result:", result)

    # question = "List protein complexes."
    # result = rag.query(question)
    # print("Query Result:", result)

    print("Welcome to the Multiomics GraphRAG system!")
    print(
        "You can ask questions about genes, transcripts, proteins, variants, structure and so on. Once finished, type 'exit' to quit."
    )
    question = ""
    while True:
        try:
            question = input("> ")
            if question.lower() != "exit":
                result = rag.query(question)
                print("----------------------------------------")
                print(result.get("answer", "No answer available."))
                print("----------------------------------------")
                print(result.get("cypher_statement", "No cypher statement available."))
                print("----------------------------------------")
                print(result.get("steps", "No steps available."))
                print("\n")
            else:
                print("Exiting the Multiomics GraphRAG system. Goodbye!")
                break
        except Exception as e:
            print(f"An error occurred: {e}")
            continue
