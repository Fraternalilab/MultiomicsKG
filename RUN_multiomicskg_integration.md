- Set environment variables and especially DATASTORE value
    ```bash
    set -a ; source .env ; set +a
    ```
- Start the docker compose application `mkgint`. This will open an interactive shell with the Python environment already installed.
    ```bash
    docker compose -f docker/docker-compose-integration.yaml up -d mkgint
    docker compose -f docker/docker-compose-integration.yaml attach mkgint
    ```
- Download the sources to the datastore automatically:
    ```bash
    python3 download_sources.py -d /root/datastore
    ```
- Install manually the source files that are not public
- Install VEP
    ```bash
    cd /root/MultiomicsKG/deps/ensembl-vep
    perl INSTALL.pl -a acp -s homo_sapiens -y GRCh38 -g all --NO_HTSLIB --NO_UPDATE --NO_TEST --CACHEDIR /root/datastore/source/vep_data_116
    ```
- Run the structural data generator process:
    ```bash
    python3 structural_data_generator.py \
                -d /root/datastore \
                -s /root/datastore/scrap \
                -p 202609A_
    ```
- Run the structural data compiler process:
    ```bash
    python3 structural_data_compiler.py \
                -d /root/datastore \
                -s /root/datastore/scrap \
                -p 202609A
    ```
- Run the integration process to generate the final CSV files. Tip: change the version prefix to a more suitable one if needed.
    ```bash
    python3 multiomicskg.py \
                -d /root/datastore \
                -v /root/MultiomicsKG/deps/ensembl-vep/vep \
                -p 202609A
    ```