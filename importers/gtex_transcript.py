from parsers.icsv import ICSV
import pandas
import shutil
import scipy


class GTexTranscript:
    def __init__(
        self,
        gtex_sample_attributes_filename: str,
        gtex_transcript_tpms_filename: str,
        node_transcript_filename: str,
    ):
        self.gtex_sample_attributes_filename = gtex_sample_attributes_filename
        self.gtex_transcript_tpms_filename = gtex_transcript_tpms_filename
        self.node_transcript_filename = node_transcript_filename
        self.gtex_sample_attributes_by_tissue_type_db = dict()
        self.transcripts_df = pandas.DataFrame()

    def set_output_files(self):
        # make a back up copy first of original node_transcript file
        shutil.copy2(
            self.node_transcript_filename, self.node_transcript_filename + ".orig"
        )

        # no need to reopen output csv file - using pandas dumping instead
        # self.node_transcript_file = open(self.node_transcript_filename, "w")
        # self.node_transcript_csv = csv.writer(
        #     self.node_transcript_file, delimiter=";", quotechar="'"
        # )

    # def set_node_transcript_header(self):
    #     header_columns = ["id:ID(transcript-id)"]
    #     for tissue in self.gtex_sample_attributes_by_tissue_type_db:
    #         tissue_name = "gtex_" + tissue.lower()
    #         tissue_name = tissue_name.replace(" ", "_")
    #         header_columns.append(tissue_name + "_median")
    #         header_columns.append(tissue_name + "_mad")

    #     # self.transcripts_df = pandas.DataFrame(
    #     #     data=None, index=None, columns=header_columns
    #     # )
    #     # self.node_transcript_csv.writerow(header_columns)

    def importfile(self):
        gtex_sample_attributes_csv = ICSV(
            self.gtex_sample_attributes_filename,
            "gtex_sample_attributes",
            max_items=-1,
            delimiter="\t",
        )

        # gtex_transcript_tpms_csv = ICSV(
        #     self.gtex_transcript_tpms_filename,
        #     "gtex_transcript_tpms",
        #     max_items=-1,
        #     delimiter="\t",
        #     estimated_items=199327,
        #     compressed=True,
        #     skip_rows=2,
        # )

        # get list of samples by tissue type
        # self.gtex_sample_attributes_by_tissue_type_db = (
        #     gtex_sample_attributes_csv.get_dict_items("SMTSD")
        # )
        while True:
            item = gtex_sample_attributes_csv.get_next_item()
            if item is None:
                break

            tissue = item["SMTSD"]
            sample_id = item["SAMPID"]
            sample_included_freeze = item["SMAFRZE"]
            # What's the criteria to select samples?
            # rna_quality > 6.0 does not select all samples considered in GTEx v8 website

            # rna_quality = 0.0

            # try:
            #     rna_quality = float(item["SMRIN"])
            # except Exception:
            #     pass

            # if rna_quality >= 6.0:
            if sample_included_freeze == "RNASEQ":
                if tissue in self.gtex_sample_attributes_by_tissue_type_db:
                    self.gtex_sample_attributes_by_tissue_type_db[tissue].append(
                        sample_id
                    )
                else:
                    self.gtex_sample_attributes_by_tissue_type_db[tissue] = list()
                    self.gtex_sample_attributes_by_tissue_type_db[tissue].append(
                        sample_id
                    )

        gtex_sample_attributes_csv.close()

        # get list of selected transcripts
        node_transcript_csv = ICSV(
            self.node_transcript_filename,
            "node_transcript",
            max_items=-1,
            delimiter=";",
        )

        selected_transcripts = dict()
        while True:
            item = node_transcript_csv.get_next_item()
            if item is None:
                break

            transcript_id = item["id:ID(transcript-id)"]
            unv_transcript_id = transcript_id.split(".")[0]

            if unv_transcript_id in selected_transcripts:
                raise Exception(
                    "Multiple versions of the same basic transcript detected"
                )
            else:
                selected_transcripts[unv_transcript_id] = list()
                selected_transcripts[unv_transcript_id].append(transcript_id)

        node_transcript_csv.close()
        self.set_output_files()

        # self.set_node_transcript_header()

        # read tpm data for each transcript
        with pandas.read_csv(
            self.gtex_transcript_tpms_filename, sep="\t", skiprows=2, iterator=True
        ) as transcript_tpms_it:
            stop = False
            transcripts_read = 0
            transcripts_chunk_df = pandas.DataFrame
            while True:
                try:
                    transcripts_chunk_df = transcript_tpms_it.get_chunk(5000)

                except StopIteration:
                    stop = True

                transcripts_read = transcripts_read + len(transcripts_chunk_df.index)
                print(f"Read {transcripts_read} transcripts data samples", end="\r")
                for row_pair in transcripts_chunk_df.iterrows():
                    transcript_data_ss = row_pair[1]
                    unv_transcript_id = transcript_data_ss.transcript_id.split(".")[0]

                    if unv_transcript_id in selected_transcripts:
                        full_transcript_id = selected_transcripts[unv_transcript_id][0]
                        self.generate_transcript_record(
                            unv_transcript_id, full_transcript_id, transcript_data_ss
                        )
                        # remove augmented transcript from list of selected ones
                        selected_transcripts.pop(unv_transcript_id)
                if stop:
                    break

        # generate null data for the remaining selected transcripts without TPM data
        for (
            remaining_transcript_unv_id,
            remaining_transcript_full_id,
        ) in selected_transcripts.items():
            self.generate_transcript_record(
                remaining_transcript_unv_id, remaining_transcript_full_id, None
            )

        # save node_transcript with tpm transcripts per tissue data
        self.transcripts_df.to_csv(self.node_transcript_filename, sep=";", index=False)

        # for transcript in transcripts:
        #     print(transcript)
        #     for tissue, samples_att in gtex_sample_attributes_by_tissue_type_db:
        #         print(tissue)
        #         for sample in samples_att:
        #             print(sample["SAMPID"])

        # get list of transcripts by tissue type and set unversioned ids

        # gtex_sample_attributes_by_tissue_type_db = dict()
        # while True:
        #     item = gtex_sample_attributes_csv.get_next_item()
        #     if item is None:
        #         break

        #     tissue = item["SMTSD"]

        #     if tissue in gtex_sample_attributes_by_tissue_type_db:
        #         gtex_sample_attributes_by_tissue_type_db[tissue].append(item)
        #     else:
        #         gtex_sample_attributes_by_tissue_type_db[tissue] = list()
        #         gtex_sample_attributes_by_tissue_type_db[tissue].append(item)

    def generate_transcript_record(
        self,
        unv_transcript_id: str,
        full_transcript_id: str,
        transcript_data_ss: pandas.Series,
    ) -> dict:
        # generate transcript statistics per tissue type
        entry = dict()
        entry["id:ID(transcript-id)"] = full_transcript_id

        for tissue, samples in self.gtex_sample_attributes_by_tissue_type_db.items():
            tissue_name = "gtex_" + tissue.lower()
            tissue_name = tissue_name.replace(" ", "_")
            # header_columns.append(tissue_name + "_median")
            # header_columns.append(tissue_name + "_mad")
            if transcript_data_ss is None:
                # output null data as TPM data is not available
                entry[tissue_name + "_median:double"] = ""
                entry[tissue_name + "_mad:double"] = ""
            else:
                tissue_samples_values = transcript_data_ss.get(samples)
                if tissue_samples_values is None:
                    raise Exception("No values retrieved for transcript by tissue type")
                entry[tissue_name + "_median:double"] = tissue_samples_values.median()
                entry[tissue_name + "_mad:double"] = scipy.stats.median_abs_deviation(
                    tissue_samples_values, scale=1
                )
            # Density Plot of samples by tissue
            # tissue_samples_values.plot.kde(title=tissue)
            # print("Pause")

        # self.transcripts_df = self.transcripts_df.append(entry, ignore_index=True)
        entry_df = pandas.DataFrame(entry, index=[full_transcript_id])
        self.transcripts_df = pandas.concat(
            [self.transcripts_df, entry_df], ignore_index=False, axis=0
        )
