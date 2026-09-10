### MultiomicsKG LLM-Enabled Querying Interfaces

This file contains instructions on how to import the MultiomicsKG database into a Neo4J instance and to run the GraphRAG Chatbot agent and the MCP server for querying the KG using the assistance of LLMs.

## Requirements

This component requires the data files of the MultiomicsKG Assembly component in one of these formats:

1. CSV Files
2. Neo4J Dump
3. An unix-like terminal (either Linux or MacOS, others may also work)
4. A local installation of [Docker and Docker Compose](https://docs.docker.com/get-started/get-docker/) active and running.
5. A running local instance of [Ollama](https://ollama.com/) (required only if running GraphRAG Chatbot).

Also, this component requires a folder structure where to store Neo4J files (data, plugins, logs) and OpenWeb UI runtime files. See the section below for instructions.

## Installation

### 1. Configure and set the environment variables

The environment where the component runs requires a lot of configuration and setting environment variables are useful to simplify the commands. This includes the necessary configuration to run the Docker application using Docker Compose.
1. Open your Terminal and change to the folder where you cloned the MultiomicsKG repository.
2. Configure your local copy of the `.env` file accordingly.
    1. Set at least the following variable to your local configuration:
        ```
        DATASTORE
        ```
        Note: `DATASTORE` should point to a folder where you'd like to store runtime files (e.g., Neo4J database and OpenWebUI files)
    2. If you want to run the GraphRAG Chatbot, you need at least `OPENAI_API_KEY` set and optionally `ANTHROPIC_API_KEY`.

3. Export the variables in the `.env` file to your running Terminal session:

    ```bash
    set -a ; source .env ; set +a
    ```

> [!IMPORTANT]
> All commands from this point on require the environment varibales from the `.env` file to be exported in your terminal session.


### 1. Create runtime folders

Choose a root folder where to store the runtime files from Neo4J and OpenWeb UI. We will refer this root folder as `DATASTORE`.

Create the following folder structure with the commands further below:
```
<DATASTORE>/
├── neo4j/
│   ├── data/
│   ├── plugins/
│   └── logs/
└── openwebui/
```

```bash
mkdir -p $NEO4J_DATA
mkdir -p $NEO4J_PLUGINS
mkdir -p $NEO4J_LOGS
mkdir -p $OPEN_WEBUI
```

### 2. Import the MultiomicsKG generated files into a Neo4J Database

Here you have two options: either compile a Neo4J database from the CSV files or from a Neo4J database dump.

#### 2.1 Using CSV files

> [!NOTE]
> Incomplete section

#### 2.2 Using a Neo4J dump

You can download a MultiomicsKG Neo4J database dump file from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23752092.svg)](https://doi.org/10.5281/zenodo.23752092).


1. **Copy the database dump**:

    Download the dump file into `$NEO4J_DATA` folder.

2. **Load the database into Neo4J**:

    With the dump file in `$NEO4J_DATA`, run the following command:

    ```bash
    docker run --interactive --tty --rm -v $NEO4J_DATA:/data:rw neo4j/neo4j-admin:5.26.31 neo4j-admin database load neo4j --from-path=/data
    ```
    This command should load the MultiomicsKG dump into a fully functional Neo4J database named `neo4j` in your `$NEO4J_DATA` folder.

3. **Install required Neo4J plugins**:

    The LLM-enabled tools require the Neo4J APOC and GDS plugins. The default Neo4J docker container looks for the plugins in `$NEO4J_PLUGINS`. Install them as follows:

    ```bash
    wget https://github.com/neo4j/apoc/releases/download/5.26.32/apoc-5.26.32-core.jar -O $NEO4J_PLUGINS/apoc-5.26.32-core.jar
    wget https://github.com/neo4j/graph-data-science/releases/download/2.13.2/neo4j-graph-data-science-2.13.2.jar -O $NEO4J_PLUGINS/neo4j-graph-data-science-2.13.2.jar
    ```
### 3. **Set up local Ollama models**

The LLM-tools in this component supports a local (or remote) instance of [Ollama](https://ollama.com/) for running open-source local LLM models.

>[!IMPORTANT]
> This step is only required if you want to run the GraphRAG Chatbot. It is not required for the MCP Server so you can skip this section if you only want to run the MCP Server.

1. **Install Ollama**:

    In case you want to run Ollama on your local host, proceed to install it according to instructions [here](https://ollama.com/download).

2. **Start Ollama**:

    ```bash
    ollama serve &
    ```

3. **Pull the required LLM-models**:

    ```bash
    ollama pull $OLLAMA_MODEL
    ollama pull $OLLAMA_EMBEDDING_MODEL
    ```

### 4. **Run the querying interfaces from Docker Compose**

The following command should build the required docker images.

```bash
docker compose -f docker/docker-compose.yaml build
```

#### Running only the MCP Server

In this scenario, the application does not require Ollama or OpenAI/Anthropic API keys set in the environment. Start the application with this command:

```
docker compose -f docker/docker-compose.yaml up neo4j-mcp
```

#### Running only the GraphRAG Chatbot

In this scenario, the application requires Ollama running and at least an OpenAI API key set in the environment. Start the application with this command:

```
docker compose -f docker/docker-compose.yaml up open-webui
```

#### Running the full application

If you want to run the full application (GraphRAG Chatbot + MCP Server), you must to have Ollama ready and at least an OpenAI API key set in the environment. Run the following command:

```bash
docker compose -f docker/docker-compose.yaml up
```
>[!TIP]
> You may want to type "D" to detach the screen and return to the terminal. Or you can open another terminal, but remember to reload the variables from the `.env` file in the repository.

>[!NOTE]
> It may take up to 10 minutes to run the full stack. The `pipelines` container which is part of the GraphRAG Chatbot application takes a lot of time to load so it is normal if it appears to hang. You can monitor the logs of the containers with the following command if needed:
>```bash
>docker compose -f docker/docker-compose.yaml logs -f
>```

> [!IMPORTANT]
> Once the stack is running, the OpenWeb UI server, MCP server and Neo4J Browser are available in the following URLs:
>
> **OpenWebUI**: `http://localhost:3000/`
>
> **MCP Server**: `http://localhost:8000/`
>
> **Neo4J Browser**: `http://localhost:7474/`
>

### 5. **Configure OpenWebUI to access GraphRAG Chatbot**

1. Go to OpenWeb UI at `http://localhost:3000/` on your Web Browser and create your admin user.
2. Connect the GraphRAG Chatbot pipeline to Open WebUI (as an admin user):
    - Navigate to the `Settings > Admin Settings > Connections > OpenAI` API section in Open WebUI.
    - In "Manage OpenAI API Connections", click "+" to Add Connection and set the base API URL to `http://pipelines:9099` and the API key to `0p3n-w3bu!`. Save your changes. The pipelines that implement the GraphRAG Chatbot should now be active and available on the model list in the main page.

3. Manage configurations (optional):
    - In the admin panel of OpenWeb UI, go to `Admin Settings > Pipelines` tab.
    - Select your desired pipeline (e.g., `pipeline_chatgpt_rag`) and modify the valve values directly from OpenWebUI. For example, you can change the provider's model name on the fly.

4. You should be ready now to submit your prompts in a session from the main page of OpenWebUI.

### 6. **Accessing the MCP Server**

The MultiomicsKG MCP Server should be available from `http://localhost:8000/`. You should be able to configure an MCP connector to the MultiomicsKG MCP Server following the instructions of your MCP-compatible application and using these default parameters:

- Host: `http://localhost:8000/`
- Full URL: `http://localhost:8000/mcp/`
- Transport: `http`
- Authorisation Header: `"Authorization: Basic $NEO4J_BASIC_AUTH"`

#### Example: configuring Claude Desktop/Claude Code

##### Requirements

1. Claude Desktop or Claude Code installed
2. Claude CLI tool
    2.1 On MacOS, install with `brew install --cask claude-code` if you use Homebrew, or go [here](https://code.claude.com/docs/en/quickstart#step-1-install-claude-code) for manual installation.

##### Instructions

1. Create a folder for your Claude Desktop/Claude Code project
2. Test if the MCP Server running locally is accessible:
    ```bash
    curl -i http://localhost:8000/mcp
    ```
    You should get something like `Method Not Allowed: only POST is supported on /mcp`. This is expected.
3. Add the MCP server to your Claude's configuration:
    ```bash
    claude mcp add --transport http multiomicskg http://localhost:8000/mcp/ --header "Authorization: Basic $NEO4J_BASIC_AUTH"
    ```
4. Verify if the MCP server is properly configured and ready:
    ```bash
    claude mcp list
    ```
    You should check that the last line looks like:
    `multiomicskg: http://localhost:8000/mcp/ (HTTP) - ✔ Connected`
5. Open your project folder in Claude Desktop/Claude Code
6. Verify that Claude is running in your project’s folder
7. Ask your questions with instructions like “Use the MultiomicsKG MCP server….” 

> [!IMPORTANT]
> All `claude` commands must be run in your project folder. Claude saves the MCP configuration per project rather than globally.

### 7. **Stopping the querying interfaces**

Run the following command to stop the querying interfaces using Docker Composer:

```bash
docker compose -f docker/docker-compose.yaml down
```