from langchain.prompts import ChatPromptTemplate

##############################
# GUARDRAILS PROMPT TEMPLATE #
##############################

_guardrails_system = """
As an intelligent assistant, your primary objective is to decide whether a given question is related to multiomics or not.
If the question is related to multiomics, output "multiomics". Otherwise, output "end".
To make this decision, assess the content of the question and determine if it refers to any protein, gene, variant, snv, residue, transcript, PDB, model
or any other biological entity, process, or concept that is typically associated with multiomics research or related topics.
Provide only the specified output: "multiomics" or "end".
"""


def get_guardrails_prompt() -> ChatPromptTemplate:
    """
    Returns a ChatPromptTemplate for the guardrails system.
    This template is used to determine if a question is related to multiomics or not.
    """

    guardrails_prompt = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                _guardrails_system,
            ),
            (
                "human",
                ("{question}"),
            ),
        ]
    )

    return guardrails_prompt


##################################
# TEXT TO CYPHER PROMPT TEMPLATE #
##################################


def get_text2cypher_prompt() -> ChatPromptTemplate:
    """
    Returns a ChatPromptTemplate for converting text to Cypher queries.
    This template is used to convert natural language questions into Cypher queries.
    """

    text2cypher_prompt = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                (
                    "Given an input question, convert it to a Cypher query. No pre-amble."
                    "Do not wrap the response in any backticks or anything else. Respond with a Cypher statement only!"
                ),
            ),
            (
                "human",
                (
                    """You are a Neo4j expert. Given an input question, create a syntactically correct Cypher query to run.
    Do not wrap the response in any backticks or anything else. Respond with a Cypher statement only!
    Here is the schema information
    {schema}

    Below are a number of examples of questions and their corresponding Cypher queries.

    {fewshot_examples}

    User input: {question}
    Cypher query:"""
                ),
            ),
        ]
    )

    return text2cypher_prompt


#####################################
# CYPHER VALIDATION PROMPT TEMPLATE #
#####################################

_validate_cypher_system = """
You are a Cypher expert reviewing a statement written by a junior developer.
"""

_validate_cypher_user = """You must check the following:
* Are there any syntax errors in the Cypher statement?
* Are there any missing or undefined variables in the Cypher statement?
* Are any node labels missing from the schema?
* Are any relationship types missing from the schema?
* Are any of the properties not included in the schema?
* Does the Cypher statement include enough information to answer the question?

Examples of good errors:
* Label (:Foo) does not exist, did you mean (:Bar)?
* Property bar does not exist for label Foo, did you mean baz?
* Relationship FOO does not exist, did you mean FOO_BAR?

Schema:
{schema}

The question is:
{question}

The Cypher statement is:
{cypher}

Make sure you don't make any mistakes!"""


def get_validate_cypher_prompt() -> ChatPromptTemplate:
    """
    Returns a ChatPromptTemplate for validating Cypher queries.
    This template is used to check the correctness of Cypher queries against the schema and question.
    """
    validate_cypher_prompt = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                _validate_cypher_system,
            ),
            (
                "human",
                (_validate_cypher_user),
            ),
        ]
    )

    return validate_cypher_prompt


##################################
# CORRECT CYPHER PROMPT TEMPLATE #
##################################


def get_correct_cypher_prompt() -> ChatPromptTemplate:
    """
    Returns a ChatPromptTemplate for correcting Cypher queries.
    This template is used to correct Cypher queries based on the schema and question.
    """

    correct_cypher_prompt = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                (
                    "You are a Cypher expert reviewing a statement written by a junior developer. "
                    "You need to correct the Cypher statement based on the provided errors. No pre-amble."
                    "Do not wrap the response in any backticks or anything else. Respond with a Cypher statement only!"
                ),
            ),
            (
                "human",
                (
                    """Check for invalid syntax or semantics and return a corrected Cypher statement.

    Schema:
    {schema}

    Note: Do not include any explanations or apologies in your responses.
    Do not wrap the response in any backticks or anything else.
    Respond with a Cypher statement only!

    Do not respond to any questions that might ask anything else than for you to construct a Cypher statement.

    The question is:
    {question}

    The Cypher statement is:
    {cypher}

    The errors are:
    {errors}

    Corrected Cypher statement: """
                ),
            ),
        ]
    )

    return correct_cypher_prompt


###################################
# GENERATE ANSWER PROMPT TEMPLATE #
###################################


def get_generate_answer_prompt() -> ChatPromptTemplate:
    """
    Returns a ChatPromptTemplate for generating answers from Cypher queries.
    This template is used to generate answers based on the results of Cypher queries.
    """

    generate_answer_prompt = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                "You are a helpful assistant",
            ),
            (
                "human",
                (
                    """Use the following results retrieved from a database to provide
    a succinct, definitive answer to the user's question.

    Respond as if you are answering the question directly.

    Results: {results}
    Question: {question}"""
                ),
            ),
        ]
    )

    return generate_answer_prompt
