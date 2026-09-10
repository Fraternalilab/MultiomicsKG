from langchain_neo4j import Neo4jGraph
from langchain_neo4j import Neo4jVector
from langchain_neo4j.chains.graph_qa.cypher_utils import CypherQueryCorrector, Schema
from langchain_openai import OpenAIEmbeddings
from langchain_ollama import OllamaEmbeddings
from langchain_core.output_parsers import StrOutputParser
from langchain_core.example_selectors import SemanticSimilarityExampleSelector

from typing import Literal, Optional, List

from langgraph.graph import END, START, StateGraph

from neo4j.exceptions import CypherSyntaxError

from pydantic import BaseModel, Field
import chatbot_graphrag.mgr.fewshot_examples as mgrfe
import chatbot_graphrag.mgr.prompt_templates as mgrpt
from chatbot_graphrag.mgr.states import InputState, OverallState, OutputState

#################################
# QA Application with LangGraph #
#################################

# Source: https://python.langchain.com/docs/tutorials/graph/
# TODO: Reorganize the code in smaller files


def setup_langgraph(
    llm,
    enhanced_graph: Neo4jGraph,
    refresh_schema: bool,
    provider: str,
) -> StateGraph:
    """
    Setup the LangGraph workflow for MultiOmics RAG.
    """

    # Refresh schema
    if refresh_schema:
        print("Refreshing schema")
        enhanced_graph.refresh_schema()

    ###################
    # GUARDRAIL CHAIN #
    ###################

    class GuardrailsOutput(BaseModel):
        decision: Literal["multiomics", "end"] = Field(
            description="Decision on whether the question is related to multiomics"
        )

    guardrails_chain = mgrpt.get_guardrails_prompt() | llm.with_structured_output(
        GuardrailsOutput
    )

    def guardrails(state: InputState) -> OverallState:
        """
        Decides if the question is related to multiomics or not.
        """
        guardrails_output = guardrails_chain.invoke({"question": state.get("question")})
        database_records = None
        if guardrails_output.decision == "end":
            database_records = "This questions is not about multiomics or related. Therefore I cannot answer this question."
        return {
            "next_action": guardrails_output.decision,
            "database_records": database_records,
            "steps": ["guardrail"],
        }

    ######################
    # FEW-SHOT PROMPTING #
    ######################

    example_selector = None

    if provider == "openai":
        neo4j_vector_options = {
            "node_label": "fewshot_examples_openai",
            "pre_delete_collection": True,
            "index_name": "vector_fewshot_examples_openai",
        }
        example_selector = SemanticSimilarityExampleSelector.from_examples(
            mgrfe.get_fewshot_examples(),
            OpenAIEmbeddings(),
            Neo4jVector,
            k=5,
            input_keys=["question"],
            **neo4j_vector_options,
        )
    elif provider == "ollama":
        neo4j_vector_options = {
            "node_label": "fewshot_examples_ollama",
            "pre_delete_collection": True,
            "index_name": "vector_fewshot_examples_ollama",
        }
        example_selector = SemanticSimilarityExampleSelector.from_examples(
            mgrfe.get_fewshot_examples(),
            OllamaEmbeddings(model="mxbai-embed-large"),
            Neo4jVector,
            k=5,
            input_keys=["question"],
            **neo4j_vector_options,
        )
    elif provider == "anthropic":
        neo4j_vector_options = {
            "node_label": "fewshot_examples_anthropic",
            "pre_delete_collection": True,
            "index_name": "vector_fewshot_examples_anthropic",
        }
        example_selector = SemanticSimilarityExampleSelector.from_examples(
            mgrfe.get_fewshot_examples(),
            OpenAIEmbeddings(),
            Neo4jVector,
            k=5,
            input_keys=["question"],
            **neo4j_vector_options,
        )
    else:
        raise ValueError(
            f"Unsupported provider: {provider}. Supported providers are 'openai', 'ollama' and 'anthropic'."
        )

    ###########################
    # CYPHER GENERATION CHAIN #
    ###########################

    text2cypher_chain = mgrpt.get_text2cypher_prompt() | llm | StrOutputParser()

    def generate_cypher(state: OverallState) -> OverallState:
        """
        Generates a cypher statement based on the provided schema and user input
        """
        NL = "\n"
        fewshot_examples = (NL * 2).join(
            [
                f"Question: {el['question']}{NL}Cypher:{el['query']}"
                for el in example_selector.select_examples(
                    {"question": state.get("question")}
                )
            ]
        )
        generated_cypher = text2cypher_chain.invoke(
            {
                "question": state.get("question"),
                "fewshot_examples": fewshot_examples,
                "schema": enhanced_graph.schema,
            }
        )
        return {"cypher_statement": generated_cypher, "steps": ["generate_cypher"]}

    ###########################
    # CYPHER VALIDATION CHAIN #
    ###########################

    class Property(BaseModel):
        """
        Represents a filter condition based on a specific node property in a graph in a Cypher statement.
        """

        node_label: str = Field(
            description="The label of the node to which this property belongs."
        )
        property_key: str = Field(description="The key of the property being filtered.")
        property_value: str = Field(
            description="The value that the property is being matched against."
        )

    class ValidateCypherOutput(BaseModel):
        """
        Represents the validation result of a Cypher query's output,
        including any errors and applied filters.
        """

        errors: Optional[List[str]] = Field(
            description="A list of syntax or semantical errors in the Cypher statement. Always explain the discrepancy between schema and Cypher statement"
        )
        filters: Optional[List[Property]] = Field(
            description="A list of property-based filters applied in the Cypher statement."
        )

    validate_cypher_chain = (
        mgrpt.get_validate_cypher_prompt()
        | llm.with_structured_output(ValidateCypherOutput)
    )

    ##########################
    # CYPHER CORRECTOR CHAIN #
    ##########################

    # Cypher query corrector is experimental
    corrector_schema = [
        Schema(el["start"], el["type"], el["end"])
        for el in enhanced_graph.structured_schema.get("relationships")
    ]
    cypher_query_corrector = CypherQueryCorrector(corrector_schema)

    ##########################
    # CYPHER VALIDATION STEP #
    ##########################

    def validate_cypher(state: OverallState) -> OverallState:
        """
        Validates the Cypher statements and maps any property values to the database.
        """
        errors = []
        mapping_errors = []
        # Check for syntax errors
        try:
            enhanced_graph.query(f"EXPLAIN {state.get('cypher_statement')}")
        except CypherSyntaxError as e:
            errors.append(e.message)
        # Experimental feature for correcting relationship directions
        corrected_cypher = cypher_query_corrector(state.get("cypher_statement"))
        if not corrected_cypher:
            errors.append("The generated Cypher statement doesn't fit the graph schema")
        if not corrected_cypher == state.get("cypher_statement"):
            print("Relationship direction was corrected")
        # Use LLM to find additional potential errors and get the mapping for values
        llm_output = validate_cypher_chain.invoke(
            {
                "question": state.get("question"),
                "schema": enhanced_graph.schema,
                "cypher": state.get("cypher_statement"),
            }
        )
        if llm_output.errors:
            errors.extend(llm_output.errors)
        # if llm_output.filters:
        #     for filter in llm_output.filters:
        #         # Do mapping only for string values
        #         if (
        #             not [
        #                 prop
        #                 for prop in enhanced_graph.structured_schema["node_props"][
        #                     filter.node_label
        #                 ]
        #                 if prop["property"] == filter.property_key
        #             ][0]["type"]
        #             == "STRING"
        #         ):
        #             continue
        #         mapping = enhanced_graph.query(
        #             f"MATCH (n:{filter.node_label}) WHERE toLower(n.`{filter.property_key}`) = toLower($value) RETURN 'yes' LIMIT 1",
        #             {"value": filter.property_value},
        #         )
        #         if not mapping:
        #             print(
        #                 f"Missing value mapping for {filter.node_label} on property {filter.property_key} with value {filter.property_value}"
        #             )
        #             mapping_errors.append(
        #                 f"Missing value mapping for {filter.node_label} on property {filter.property_key} with value {filter.property_value}"
        #             )
        if mapping_errors:
            next_action = "end"
        elif errors:
            next_action = "correct_cypher"
        else:
            next_action = "execute_cypher"

        return {
            "next_action": next_action,
            "cypher_statement": corrected_cypher,
            "cypher_errors": errors,
            "steps": ["validate_cypher"],
        }

    correct_cypher_chain = mgrpt.get_correct_cypher_prompt() | llm | StrOutputParser()

    def correct_cypher(state: OverallState) -> OverallState:
        """
        Correct the Cypher statement based on the provided errors.
        """
        corrected_cypher = correct_cypher_chain.invoke(
            {
                "question": state.get("question"),
                "errors": state.get("cypher_errors"),
                "cypher": state.get("cypher_statement"),
                "schema": enhanced_graph.schema,
            }
        )

        return {
            "next_action": "validate_cypher",
            "cypher_statement": corrected_cypher,
            "steps": ["correct_cypher"],
        }

    ##########################
    # CYPHER EXECUTION CHAIN #
    ##########################

    no_results = "I couldn't find any relevant information in the database"

    def execute_cypher(state: OverallState) -> OverallState:
        """
        Executes the given Cypher statement.
        """

        records = enhanced_graph.query(state.get("cypher_statement"))
        return {
            "database_records": records if records else no_results,
            "next_action": "end",
            "steps": ["execute_cypher"],
        }

    #########################
    # GENERATE ANSWER CHAIN #
    #########################

    generate_final_chain = mgrpt.get_generate_answer_prompt() | llm | StrOutputParser()

    def generate_final_answer(state: OverallState) -> OutputState:
        """
        Decides if the question is related to movies.
        """
        final_answer = generate_final_chain.invoke(
            {
                "question": state.get("question"),
                "results": state.get("database_records"),
            }
        )
        return {"answer": final_answer, "steps": ["generate_final_answer"]}

    ######################
    # LANGGRAPH WORKFLOW #
    ######################

    def guardrails_condition(
        state: OverallState,
    ) -> Literal["generate_cypher", "generate_final_answer"]:
        if state.get("next_action") == "end":
            return "generate_final_answer"
        elif state.get("next_action") == "multiomics":
            return "generate_cypher"

    def validate_cypher_condition(
        state: OverallState,
    ) -> Literal["generate_final_answer", "correct_cypher", "execute_cypher"]:
        if state.get("next_action") == "end":
            return "generate_final_answer"
        elif state.get("next_action") == "correct_cypher":
            return "correct_cypher"
        elif state.get("next_action") == "execute_cypher":
            return "execute_cypher"

    ###################
    # BUILD LANGGRAPH #
    ###################

    langgraph = StateGraph(OverallState, input=InputState, output=OutputState)
    langgraph.add_node(guardrails)
    langgraph.add_node(generate_cypher)
    langgraph.add_node(validate_cypher)
    langgraph.add_node(correct_cypher)
    langgraph.add_node(execute_cypher)
    langgraph.add_node(generate_final_answer)

    langgraph.add_edge(START, "guardrails")
    langgraph.add_conditional_edges(
        "guardrails",
        guardrails_condition,
    )
    langgraph.add_edge("generate_cypher", "validate_cypher")
    langgraph.add_conditional_edges(
        "validate_cypher",
        validate_cypher_condition,
    )
    langgraph.add_edge("execute_cypher", "generate_final_answer")
    langgraph.add_edge("correct_cypher", "validate_cypher")
    langgraph.add_edge("generate_final_answer", END)

    langgraph = langgraph.compile()
    # View
    # display(Image(self.langgraph.get_graph().draw_mermaid_png()))

    return langgraph
