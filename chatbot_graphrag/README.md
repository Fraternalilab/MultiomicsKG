

poetry install
poetry env use python3.12
eval $(poetry env activate)

export PIPELINES_DIR="/Users/julian/dev/projects/phd/snvpred/openwebui_integration/pipelines/src/pipelines/"

export $(xargs <.env)

cd pipelines/upstream
./start.sh


cd openwebui_integration/openwebui
eval $(poetry env activate)
open-webui serve

Settings > Connections > OpenAI API
Set the API URL to http://localhost:9099 and the API key to 0p3n-w3bu!
In the admin panel, go to Admin Settings > Pipelines tab.
Upload pipeline



Install APOC and GDS