"""Cheap reproducible BLAST retention diagnostic, not a production rerun.

python scripts/sanity_short_queries.py OUTPUT --bases N --sequences N
"""
import argparse
import json
import random
import shutil
import subprocess
from pathlib import Path

from database_residue_screen.database import inspect_database
from database_residue_screen.example import save_fasta, save_json
from database_residue_screen.length_aware import engine_statistics, score_calibration
from database_residue_screen.matching import BLAST_FIELDS, strict_hsp, tolerant_hsp
from database_residue_screen.search import iter_hsps, search_settings
from database_residue_screen.sequences import accession_version, reverse_complement


def run(destination, bases, sequences, random_bases=3000000, diagnostic_raised_thresholds=False):
    destination.mkdir(parents=True, exist_ok=False)
    rng = random.Random(19019)
    lengths = [8, 10, 16, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 32]
    queries = {f'q{length}': ''.join(rng.choices('ACGT', k=length)) for length in lengths}
    subjects = {f'RANDOM{i:04d}.1': ''.join(rng.choices('ACGT', k=random_bases // 20)) for i in range(20)}
    truth = []
    for name, seq in queries.items():
        length = len(seq)
        # All permitted spans, each maximal error count, with internal errors.
        variants = [('exact', seq), ('reverse', reverse_complement(seq))]
        for span in range((9 * length + 9)//10, length + 1):
            clipped = seq[:span]
            variants.append((f'span{span}_exact', clipped))
            errors = span // 20
            if errors:
                mutated = list(clipped)
                for j in range(errors):
                    position = (j + 1) * span // (errors + 1)
                    mutated[position] = next(base for base in 'ACGT' if base != mutated[position])
                variants.append((f'span{span}_mismatch{errors}', ''.join(mutated)))
        for index, (variant, value) in enumerate(variants):
            identifier = f'PLANT{length:03d}{index:03d}.1'
            subjects[identifier] = 'N'*50 + value + 'N'*50
            truth.append(dict(query=name, subject=identifier, variant=variant, strict=variant in {'exact','reverse'} or value == seq,
                              tolerant=True))
    save_fasta(destination/'queries.fasta', queries)
    save_fasta(destination/'subjects.fasta', subjects)
    binary = shutil.which('blastn')
    subprocess.run(['makeblastdb','-in',str(destination/'subjects.fasta'),'-out',str(destination/'db'),'-dbtype','nucl','-parse_seqids'],
                   check=True, capture_output=True)
    actual = inspect_database(destination/'db', destination/'database_evidence', hash_all=True)
    stats, witnesses = engine_statistics(binary, lengths)
    summary, calibration = [], []
    for length in lengths:
        query = {f'q{length}': queries[f'q{length}']}
        qpath = destination/f'q{length}.fasta'
        save_fasta(qpath, query)
        calibrated = score_calibration(length, bases, sequences, stats)
        calibration.append(calibrated)
        for arm in (('production_E10', 'diagnostic_chance_saturated') if diagnostic_raised_thresholds else ('production_E10',)):
            diagnostic = score_calibration(length, bases, sequences, stats, mode='adaptive')
            settings = search_settings(2000, 10 if arm == 'production_E10' else diagnostic['evalue'], calibrated['effective_search_space'])
            output = destination/f'{arm}_{length}.tsv'
            command = [binary,'-db',str(destination/'db'),'-query',str(qpath),'-outfmt','7 '+' '.join(BLAST_FIELDS),'-out',str(output),'-num_threads','1']
            for key,value in settings.items():
                if value is not None:
                    command += ['-'+key,str(value)]
            result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=60)
            hsps = list(iter_hsps(output, query))
            # accession_version handles lcl and accession wrappers consistently.
            observed = {accession_version(h['sseqid']) for h in hsps if tolerant_hsp(h)}
            expected = [row for row in truth if row['query'] in query]
            summary.append(dict(length=length, arm=arm, strict_testability=calibrated['strict_testability'],
                                tolerant_testability=calibrated['tolerant_testability'], partially_score_censored=calibrated['partially_score_censored'], evalue=settings['evalue'], effective_search_space=settings['searchsp'],
                                command=command, stderr=result.stderr, planted_expected=len(expected),
                                planted_tolerant_recovered=sum(row['subject'] in observed for row in expected),
                                planted_exact_expected=sum(row['strict'] for row in expected),
                                planted_exact_recovered=sum(row['strict'] and row['subject'] in observed for row in expected),
                                raw_hsps=len(hsps), random_hsps=sum('RANDOM' in h['sseqid'] for h in hsps),
                                random_strict_hsps=sum('RANDOM' in h['sseqid'] and strict_hsp(h) for h in hsps),
                                random_tolerant_hsps=sum('RANDOM' in h['sseqid'] and tolerant_hsp(h) for h in hsps),
                                missing_plants=[row['variant'] for row in expected if row['subject'] not in observed]))
    report = dict(seed=19019, blast_version=subprocess.check_output([binary,'-version'],text=True).strip(),
                  actual_database=actual['runtime_metadata'], emulated_database=dict(bases=bases,sequences=sequences),
                  random_background_bases=random_bases, calibration=calibration, local_engine_witnesses=witnesses,
                  truth=truth, results=summary, scope='Synthetic diagnostic with per-length -searchsp; no core_nt search and no production recall claim')
    save_json(destination/'sanity.json',report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination',type=Path)
    parser.add_argument('--bases',type=int,required=True)
    parser.add_argument('--sequences',type=int,required=True)
    parser.add_argument('--diagnostic-raised-thresholds', action='store_true', help='Optional chance-saturated diagnostic, default OFF')
    args=parser.parse_args()
    result=run(args.destination.resolve(),args.bases,args.sequences,diagnostic_raised_thresholds=args.diagnostic_raised_thresholds)
    print(json.dumps(result['results'],indent=2))
