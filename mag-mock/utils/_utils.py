import glob
import sys
import os
import numpy as np
import pandas as pd
import skbio
import subprocess
import tempfile

import qiime2 as q2
from qiime2.plugins import (
    quality_control,
)



def generate_abundances(profile, num_genomes, mu=0, sigma=1, lambd=0.5):
    abundances = []
    if profile == "uniform":
        abundances = [1/num_genomes] * num_genomes
    elif profile == "exponential":
        abundances = np.exp(-lambd * np.arange(num_genomes))
        abundances /= abundances.sum()
    elif profile == "lognormal":
        abundances = np.exp(mu + sigma * np.random.randn(num_genomes))
        abundances /= abundances.sum()
    else:
        print(f"Invalid abundance profile option: {profile}")
        sys.exit(1)
    return abundances


def abundances_to_df(abundances, genome_files, sample_id):
    data = {
        "id": [os.path.basename(f).split(".")[0] for f in genome_files],
        sample_id: abundances
    }
    df = pd.DataFrame(data)
    df.set_index("id", inplace=True)
    return df

def generate_reads(genome_files, abundances, total_reads, sample_name, tmp_dir, threads=1, read_len=150, seed=0):
    for genome_file, abundance in zip(genome_files, abundances):
        genome_reads = int(total_reads * abundance)
        _id = os.path.basename(genome_file).replace('.fasta', '')
        print(f"Generating {genome_reads} reads for {_id}")

        cmd = [
            "mason_simulator",
            "-v",
            "--seed", str(seed),
            "--illumina-read-length", str(read_len),
            "--seq-technology", "illumina",
            "--num-threads", str(threads),
            "-ir", genome_file,
            "-n", str(genome_reads),
            "-o", os.path.join(tmp_dir, f"{sample_name}_{_id}_L001_R1_001.fastq.gz"),
            "-or", os.path.join(tmp_dir, f"{sample_name}_{_id}_L001_R2_001.fastq.gz")
        ]

        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"Error generating reads for {_id}")

def combine_reads(sample_name, tmp_dir, dst_dir):
    forward_reads = os.path.join(dst_dir, f"{sample_name}_00_L001_R1_001.fastq.gz")
    reverse_reads = os.path.join(dst_dir, f"{sample_name}_00_L001_R2_001.fastq.gz")
    
    with open(forward_reads, 'wb') as f_out:
        for file in sorted(os.listdir(tmp_dir)):
            if file.endswith("_L001_R1_001.fastq.gz"):
                with open(os.path.join(tmp_dir, file), 'rb') as f_in:
                    f_out.write(f_in.read())
    
    with open(reverse_reads, 'wb') as f_out:
        for file in sorted(os.listdir(tmp_dir)):
            if file.endswith("_L001_R2_001.fastq.gz"):
                with open(os.path.join(tmp_dir, file), 'rb') as f_in:
                    f_out.write(f_in.read())

    print(f"Combining all forward reads into {sample_name}_00_L001_R1_001.fastq.gz")
    print(f"Combining all reverse reads into {sample_name}_00_L001_R2_001.fastq.gz")

def simulate_reads(
    genomes_dir, total_reads, abundance_profile, sample_name, 
    mu=0, sigma=1, lambd=0.5, simulated_reads_dir='.', threads=1, 
    read_len=150, seed=0
):
    genome_files = glob.glob(os.path.join(genomes_dir, "*.fasta"))
    num_genomes = len(genome_files)

    if num_genomes == 0:
        print("No genome files found in the specified directory.")
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmp:
        abundances = generate_abundances(abundance_profile, num_genomes, mu, sigma, lambd)
        abundances_df = abundances_to_df(abundances, genome_files, sample_name)
        generate_reads(
            genome_files=genome_files, abundances=abundances, total_reads=total_reads, 
            sample_name=sample_name, tmp_dir=tmp, threads=threads, read_len=read_len, 
            seed=seed
        )
        combine_reads(sample_name, tmp, simulated_reads_dir)
        print(os.listdir(simulated_reads_dir))

    print(f"Simulated reads are saved as {sample_name}_00_L001_R1_001.fastq.gz and {sample_name}_00_L001_R2_001.fastq.gz")
    return abundances_df

def krakenize_taxonomy(taxonomy: pd.Series) -> pd.Series:
    taxa = taxonomy.copy()
    
    # add domain information - careful! only valid for Bacteria
    taxa = "d__Bacteria;" + taxa

    # insert genus information into species name
    def add_genus(t: str) -> str:
        t = t.split(";")
        genus = t[-2].split("__")[-1]
        species = t[-1].split("__")[-1]
        new_species = f"s__{genus} {species}"
        new_t = ";".join(t[:-2])
        new_t += f";{new_species}"
        return new_t
        
    taxa = taxa.apply(lambda x: add_genus(x))
    return taxa


def get_human_ids(fp_fwd, fp_rev):
    ids_fwd = [
        r.metadata["id"] for r in skbio.io.read(fp_fwd, format="fastq", phred_offset=33)
    ]
    ids_rev = [
        r.metadata["id"] for r in skbio.io.read(fp_rev, format="fastq", phred_offset=33)
    ]
    return ids_fwd, ids_rev


def spike_reads(sample, sample_name, human_fwd, human_rev):
    with tempfile.TemporaryDirectory() as tmp:
        sample.export_data(tmp)
        fwd_fname = f"{sample_name}_00_L001_R1_001.fastq.gz"
        fwd = os.path.join(tmp, fwd_fname)
        rev_fname = f"{sample_name}_00_L001_R2_001.fastq.gz"
        rev = os.path.join(tmp, rev_fname)
    
        samples_new_dir = os.path.join(tmp, "samples_new")
        os.makedirs(samples_new_dir, exist_ok=True)
        
        sf_new_fp = os.path.join(samples_new_dir, fwd_fname)
        with open(sf_new_fp, "ab") as fout:
            p1 = subprocess.Popen(["cat", fwd], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["cat", human_fwd], stdout=subprocess.PIPE)
            fout.write(p1.communicate()[0])
            fout.write(p2.communicate()[0])
    
        sr_new_fp = os.path.join(samples_new_dir, rev_fname)
        with open(sr_new_fp, "ab") as fout:
            p1 = subprocess.Popen(["cat", rev], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["cat", human_rev], stdout=subprocess.PIPE)
            fout.write(p1.communicate()[0])
            fout.write(p2.communicate()[0])
    
        return q2.Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]",
            samples_new_dir,
            "CasavaOneEightSingleLanePerSampleDirFmt",
        )


def _filter_reads_native(demultiplexed_sequences, database, n_threads, sensitivity, sample_name):
    with tempfile.TemporaryDirectory() as tmp:
        demultiplexed_sequences.export_data(tmp)
        fwd_fname = f"{sample_name}_00_L001_R1_001.fastq.gz"
        fwd = os.path.join(tmp, fwd_fname)
        rev_fname = f"{sample_name}_00_L001_R2_001.fastq.gz"
        rev = os.path.join(tmp, rev_fname)
        
        db_dir = os.path.join(tmp, "db")
        database.export_data(db_dir)
        print(os.listdir(db_dir))
        database_files = glob.glob(os.path.join(db_dir, "*.bt2*"))
        print("Index files:", database_files)
        index_prefix = database_files[0].split(".")[0]

        sam_out = os.path.join(tmp, "alignment.sam")
        cmd = ["bowtie2", "-p", str(n_threads), f"--{sensitivity}-local", "--rfg", "5,3",
               "-x", index_prefix, "-1", fwd, "-2", rev, "-S", sam_out]
        print(f"Running {' '.join(cmd)}")
        subprocess.run(cmd)
        
        bam_out_keep = os.path.join(tmp, "samtools_keep.bam")
        cmd = ["samtools", "view", "-b", sam_out, "-o", bam_out_keep, "-@", str(n_threads), "-F", "256", "-f", "12"]
        print(f"Running {' '.join(cmd)}")
        subprocess.run(cmd)
        
        bam_out_discard = os.path.join(tmp, "samtools_discard.bam")
        cmd = ["samtools", "view", "-b", sam_out, "-o", bam_out_discard, "-@", str(n_threads), "-F", "268"]
        print(f"Running {' '.join(cmd)}")
        subprocess.run(cmd)
        
        fastq_dir_keep = os.path.join(tmp, "fastq_keep")
        fastq_dir_discard = os.path.join(tmp, "fastq_discard")
        os.makedirs(fastq_dir_keep, exist_ok=True)
        os.makedirs(fastq_dir_discard, exist_ok=True)
        
        fastq_keep_fwd = os.path.join(fastq_dir_keep, f"{sample_name}_00_L001_R1_001.fastq.gz")
        fastq_keep_rev = os.path.join(fastq_dir_keep, f"{sample_name}_00_L001_R2_001.fastq.gz")
        with tempfile.NamedTemporaryFile() as out:
            cmd = ["samtools", "sort", "-n", "-@", str(n_threads - 1), "-o", out.name, bam_out_keep]
            print(f"Running {' '.join(cmd)}")
            subprocess.run(cmd)
            cmd = ["samtools", "fastq", "-0", "/dev/null", "-1", fastq_keep_fwd, "-2", fastq_keep_rev, "-@", str(n_threads - 1), "-n", bam_out_keep]
            print(f"Running {' '.join(cmd)}")
            subprocess.run(cmd)
        
        fastq_discard_fwd = os.path.join(fastq_dir_discard, f"{sample_name}_00_L001_R1_001.fastq.gz")
        fastq_discard_rev = os.path.join(fastq_dir_discard, f"{sample_name}_00_L001_R2_001.fastq.gz")
        with tempfile.NamedTemporaryFile() as out:
            cmd = ["samtools", "sort", "-n", "-@", str(n_threads - 1),"-o", out.name, bam_out_discard]
            print(f"Running {' '.join(cmd)}")
            subprocess.run(cmd)
            cmd = ["samtools", "fastq", "-0", "/dev/null", "-1", fastq_discard_fwd, "-2", fastq_discard_rev, "-@", str(n_threads - 1), "-n", bam_out_discard]
            print(f"Running {' '.join(cmd)}")
            subprocess.run(cmd)
        
        not_aligned = q2.Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]",
            fastq_dir_keep,
            "CasavaOneEightSingleLanePerSampleDirFmt",
        ) # reads which aligned
        aligned = q2.Artifact.import_data(
            "SampleData[PairedEndSequencesWithQuality]",
            fastq_dir_discard,
            "CasavaOneEightSingleLanePerSampleDirFmt",
        )
        
        return not_aligned, aligned
    

def filter_reads(reads, indices, threads, exclude=False, sensitivity="very-sensitive"):
    results= []
    for index in indices:
        filtered_out, = quality_control.methods.filter_reads(
            demultiplexed_sequences=reads,
            database=index,
            n_threads=threads,
            exclude_seqs=exclude,
            sensitivity=sensitivity
        )
        results.append(filtered_out)
    return results


def filter_reads_sequential(reads, indices, threads, sample_name):
    retained_reads = reads
    combined_reads = None
    with tempfile.TemporaryDirectory() as tmp:
        for i, index in enumerate(indices):
            print(f"Filtering against {i}th index: {index}")
            filtered_out, = quality_control.methods.filter_reads(
                demultiplexed_sequences=retained_reads,
                database=index,
                n_threads=threads,
                exclude_seqs=False,
                sensitivity="very-sensitive"
            )
            retained_reads, = quality_control.methods.filter_reads(
                demultiplexed_sequences=retained_reads,
                database=index,
                n_threads=threads,
                exclude_seqs=True,
                sensitivity="very-sensitive"
            )

            # add the filtered out reads to the final collection 
            out_dir = os.path.join(str(tmp), f"filtered_{i}")
            if combined_reads is None:
                combined_reads = filtered_out
            else:
                print("Extracting filtered reads to", out_dir)
                filtered_out.export_data(out_dir)
                fwd_read = os.path.join(out_dir, f"{sample_name}_00_L001_R1_001.fastq.gz")
                rev_read = os.path.join(out_dir, f"{sample_name}_00_L001_R2_001.fastq.gz")
                print("Combining reads")
                if os.path.getsize(fwd_read) > 0 and os.path.getsize(rev_read) > 0:
                    combined_reads = spike_reads(combined_reads, sample_name, fwd_read, rev_read)
                else:
                    print("One of {fwd_read}, {rev_read} is empty - ignoring...")
    
    return combined_reads


def filter_reads_sequential_NATIVE(reads, indices, threads, sample_name):
    retained_reads = reads
    combined_reads = None
    with tempfile.TemporaryDirectory() as tmp:
        for i, index in enumerate(indices):
            print(f"Filtering against {i}th index: {index}")
            (retained_reads, filtered_out) = _filter_reads_native(
                demultiplexed_sequences=retained_reads,
                database=index,
                n_threads=threads,
                sensitivity="very-sensitive",
                sample_name=sample_name
            )

            # add the filtered out reads to the final collection 
            out_dir = os.path.join(str(tmp), f"filtered_{i}")
            if combined_reads is None:
                combined_reads = filtered_out
            else:
                print("Extracting filtered reads to", out_dir)
                filtered_out.export_data(out_dir)
                fwd_read = os.path.join(out_dir, f"{sample_name}_00_L001_R1_001.fastq.gz")
                rev_read = os.path.join(out_dir, f"{sample_name}_00_L001_R2_001.fastq.gz")
                print("Combining reads")
                combined_reads = spike_reads(combined_reads, sample_name, fwd_read, rev_read)
    
    return combined_reads


def compare_results(
    human_sample_ids: list, human_read_ids_all: list, 
    filtering_results: dict, sample_name: str
    ):
    
    results = {}
    for index_name, filtering_result in filtering_results.items():
        results[index_name] = {}
        zipped = zip(human_sample_ids, human_read_ids_all, filtering_result)
        for human_sample_id, human_read_ids, filtering_result in zipped:
            with tempfile.TemporaryDirectory() as tmp:
                filtered_dir = os.path.join(tmp, f"filtered_{human_sample_id}")
                filtering_result.export_data(filtered_dir)
                
                fastq_fp = os.path.join(filtered_dir, f"{sample_name}_00_L001_R1_001.fastq.gz")
                ids_fwd = [r.metadata["id"] for r in skbio.io.read(fastq_fp, format="fastq", phred_offset=33)]
                
                total_reads = len(ids_fwd)
                intersect_reads = len(set(human_read_ids).intersection(ids_fwd))
                intersect_frac = 100 * intersect_reads / len(human_read_ids)
                
                results[index_name][human_sample_id] = {
                    'total_reads': total_reads,
                    'intersect_reads': intersect_reads,
                    'intersect_frac': intersect_frac,
                    'non_human': total_reads - intersect_reads
                }
    
    rows = []
    for index, samples in results.items():
        for sample, metrics in samples.items():
            row = {"type": index, "sample": sample}
            row.update(metrics)
            rows.append(row)

    return pd.DataFrame(rows)
