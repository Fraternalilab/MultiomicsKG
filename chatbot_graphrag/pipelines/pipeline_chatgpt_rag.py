import os
import logging
from pydantic import BaseModel
from typing import List, Union, Generator, Iterator
from dotenv import load_dotenv
from langchain_neo4j import Neo4jGraph
from chatbot_graphrag.multiomics_graph_rag import MultiomicsGraphRAG

# References
# https://python.langchain.com/docs/integrations/graphs/neo4j_cypher/
# https://python.langchain.com/docs/tutorials/graph/
# https://python.langchain.com/api_reference/reference.html
# https://python.langchain.com/api_reference/neo4j/chains/langchain_neo4j.chains.graph_qa.cypher.GraphCypherQAChain.html
# https://docs.openwebui.com/pipelines/
# https://github.com/open-webui/pipelines/blob/main/examples/pipelines/rag/text_to_sql_pipeline.py
# https://github.com/JordanNanos/example-pipelines/blob/main/01_text_to_sql_pipeline_vLLM_llama.py


load_dotenv()

logging.basicConfig(level=logging.DEBUG)


class Pipeline:
    class Valves(BaseModel):
        OPENAI_MODEL: str

    def __init__(self):
        # Optionally, you can set the id and name of the pipeline.
        # Best practice is to not specify the id so that it can be automatically inferred from the filename, so that users can install multiple versions of the same pipeline.
        # The identifier must be unique across all pipelines.
        # The identifier must be an alphanumeric string that can include underscores or hyphens. It cannot contain spaces, special characters, slashes, or backslashes.
        # self.id = "pipeline_example"

        self.NEO4J_URI = os.getenv("NEO4J_URI", "")
        self.NEO4J_USERNAME = os.getenv("NEO4J_USERNAME", "")
        self.NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "")
        self.NEO4J_DATABASE = os.getenv("NEO4J_DATABASE", "")
        self.OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "")
        self.OPENAI_MODEL = os.getenv("OPENAI_MODEL", "")

        self.enhanced_graph = Neo4jGraph(
            url=self.NEO4J_URI,
            username=self.NEO4J_USERNAME,
            password=self.NEO4J_PASSWORD,
            database=self.NEO4J_DATABASE,
            enhanced_schema=True,
            timeout=120,
        )

        self.valves = self.Valves(
            **{
                "pipelines": ["*"],
                "OPENAI_MODEL": self.OPENAI_MODEL,
            }
        )

        # The name of the pipeline.
        self.name = f"MultiomicsKB RAG (OpenAI: {self.valves.OPENAI_MODEL})"

    async def on_startup(self):
        # This function is called when the server is started.
        print(f"on_startup:{__name__}")

        print(f"{__name__}:openai model by valve: {self.valves.OPENAI_MODEL}")

        self.name = f"MultiomicsKB RAG (OpenAI: {self.valves.OPENAI_MODEL})"

        print(f"{__name__}:refreshing schema in Neo4j graph")
        self.enhanced_graph.refresh_schema()

        self.rag = MultiomicsGraphRAG(
            enhanced_graph=self.enhanced_graph,
            refresh_schema=False,
            provider="openai",
            model=self.valves.OPENAI_MODEL,
            api_key=self.OPENAI_API_KEY,
        )

    async def on_shutdown(self):
        # This function is called when the server is stopped.
        print(f"on_shutdown:{__name__}")
        pass

    async def on_valves_updated(self):
        # This function is called when the valves are updated.
        print(f"on_valves_updated:{__name__}")
        print(f"openai model by valve: {self.valves.OPENAI_MODEL}")
        self.name = f"MultiomicsKB RAG (OpenAI: {self.valves.OPENAI_MODEL})"
        self.rag = MultiomicsGraphRAG(
            enhanced_graph=self.enhanced_graph,
            refresh_schema=False,
            provider="openai",
            model=self.valves.OPENAI_MODEL,
            api_key=self.OPENAI_API_KEY,
        )

    async def inlet(self, body: dict, user: dict) -> dict:
        # This function is called before the OpenAI API request is made. You can modify the form data before it is sent to the OpenAI API.
        print(f"inlet:{__name__}")

        # print(body)
        # print(user)

        return body

    async def outlet(self, body: dict, user: dict) -> dict:
        # This function is called after the OpenAI API response is completed. You can modify the messages after they are received from the OpenAI API.
        print(f"outlet:{__name__}")

        # print(body)
        # print(user)

        return body

    def pipe(
        self, user_message: str, model_id: str, messages: List[dict], body: dict
    ) -> Union[str, Generator, Iterator]:
        # This is where you can add your custom pipelines like RAG.
        print(f"pipe:{__name__}")

        # # If you'd like to check for title generation, you can add the following check
        # if body.get("title", False):
        #     print("Title Generation Request")

        # print(f"messages\n-----\n{messages}")
        # print(f"user_message\n-----\n{user_message}")
        # print(f"body\n-----\n{body}")

        output = ""
        try:
            result = self.rag.query(user_message)

            cyper_query = result.get(
                "cypher_statement", "No cypher statement available."
            )
            final_result = result.get("answer", "No answer available.")

            output = f"```cypher\n{cyper_query}\n```\n\n{final_result}\n"
        except Exception as e:
            output = f"Exception:\n----------\n\n{e}"

        return output
