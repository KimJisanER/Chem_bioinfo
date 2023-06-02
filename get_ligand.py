#!/usr/bin/env python

import os
from chembl_webresource_client.settings import Settings
from chembl_webresource_client.new_client import new_client
from functools import lru_cache
from tqdm import tqdm


# Settings.Instance().CACHING = True
# @lru_cache(maxsize=128)
# def download_chembl_compounds(accession,save_path,gene_name):
#
#     # Create the "ligands" directory if it does not exist
#     if not os.path.exists(os.path.join(save_path, "ligands")):
#         os.makedirs(os.path.join(save_path, "ligands"))
#
#     targets_api = new_client.target
#     compounds_api = new_client.molecule
#     bioactivities_api = new_client.activity
#
#     targets = targets_api.get(target_components__accession=accession,
#                               target_type__in=['SINGLE PROTEIN','PROTEIN FAMILY']).only("target_chembl_id")
#     targets_id = list(set([target['target_chembl_id'] for target in targets]))
#
#     if len(targets_id) != 0:
#
#         target_columns_list = ['molecule_chembl_id', 'molecule_name', 'molecule_max_phase', 'molecular_weight', 'alogp',
#                                'standard_type', 'standard_relation', 'standard_value', 'standard_units',
#                                'pchembl_value', 'assay_type', 'target_name','target_organism', 'target_type']
#
#         mol_file = ''
#         for target in targets_id:
#
#             # bioactivities = bioactivities_api.filter(target_chembl_id=target,
#             #                                          standard_type__in=["IC50", "EC50", "Ki", "Kd"],
#             #                                          assay_type__in=["B", "F"], relation__in=['=', '<', '<=']).only(
#             #     *target_columns_list)
#
#             bioactivities = bioactivities_api.filter(target_chembl_id=target,
#                                                      standard_type__in=["IC50"],
#                                                      assay_type__in=["B", "F"], relation__in=['=', '<', '<=']).only(
#                 *target_columns_list)
#
#             bioactivities_list = list(bioactivities)
#
#             compound_columns_list = ["molecule_chembl_id", 'molecule_structures', 'molecule_properties', 'mw_freebase']
#             compounds = compounds_api.filter(
#                 molecule_chembl_id__in=[x['molecule_chembl_id'] for x in bioactivities],
#                 molecule_properties__mw_freebase__lte=700, molfile__isnull=False).only(*compound_columns_list)
#
#             compounds_list = [record for record in compounds]
#
#             if len(compounds_list) != 0:
#                 for record in tqdm(compounds_list, total=len(compounds_list), desc=f'Retrieving ligands for {accession}'):
#                     if record['molecule_structures']['molfile'] is not None:
#                         molfile = record['molecule_structures']['molfile']
#                         mw_freebase = f"> <mw_freebase>\n{record['molecule_properties']['mw_freebase']}\n"
#                         alogp = f"\n> <alogp>\n{record['molecule_properties']['alogp']}\n"
#                         mol_file += molfile
#                         mol_file += mw_freebase
#                         mol_file += alogp
#                         mol_file += "\n> <bioactivities>\n"
#
#                         pChEMBL_value = 0
#                         for i in bioactivities_list:
#                             if i['molecule_chembl_id'] == record['molecule_chembl_id']:
#                                 bioactivities_str = f"assay_type: {i['assay_type']}, pchembl_value: {i['pchembl_value']}, {i['standard_type']}{i['standard_relation']}{i['standard_value']}{i['standard_units']}\n"
#                                 mol_file += bioactivities_str
#                                 if i['pchembl_value'] is not None:
#                                     if float(i['pchembl_value']) > float(pChEMBL_value):
#                                         pChEMBL_value = str(i['pchembl_value'])
#                         if pChEMBL_value == 0 or pChEMBL_value == '0':
#                             pChEMBL_value = 'None'
#                         pChEMBL_value_str = f"\n> <pChEMBL_value>\n{pChEMBL_value}\n"
#                         mol_file += pChEMBL_value_str
#
#                         mol_file += "\n$$$$\n"
#         file_name = f"{gene_name}_{accession}_ligands.sdf"
#         file_path = os.path.join(save_path, "ligands", file_name)
#         with open(file_path, 'w') as outfile:
#             outfile.write(mol_file)

# Settings.Instance().CACHING = True
# @lru_cache(maxsize=128)
# def download_chembl_compounds(accession,save_path,gene_name):
#
#     # Create the "ligands" directory if it does not exist
#     if not os.path.exists(os.path.join(save_path, "ligands")):
#         os.makedirs(os.path.join(save_path, "ligands"))
#
#     targets_api = new_client.target
#     compounds_api = new_client.molecule
#     bioactivities_api = new_client.activity
#
#     targets = targets_api.get(target_components__accession=accession,
#                               target_type__in=['SINGLE PROTEIN','PROTEIN FAMILY']).only("target_chembl_id")
#     targets_id = list(set([target['target_chembl_id'] for target in targets]))
#
#     if len(targets_id) != 0:
#
#         target_columns_list = ['molecule_chembl_id', 'molecule_name', 'molecule_max_phase', 'molecular_weight', 'alogp',
#                                'standard_type', 'standard_relation', 'standard_value', 'standard_units',
#                                'pchembl_value', 'assay_type', 'target_name','target_organism', 'target_type']
#
#         mol_file = ''
#         for target in targets_id:
#
#             bioactivities = bioactivities_api.filter(target_chembl_id=target,
#                                                      standard_type__in=["IC50", "EC50", "Ki", "Kd"],
#                                                      assay_type__in=["B", "F"], relation__in=['=', '<', '<=']).only(
#                 *target_columns_list)
#
#
#             bioactivities_list = list(bioactivities)
#
#             compound_columns_list = ["molecule_chembl_id", 'molecule_structures', 'molecule_properties', 'mw_freebase']
#             compounds = compounds_api.filter(
#                 molecule_chembl_id__in=[x['molecule_chembl_id'] for x in bioactivities],
#                 molecule_properties__mw_freebase__lte=700, molfile__isnull=False).only(*compound_columns_list)
#
#             compounds_list = [record for record in compounds]
#
#             if len(compounds_list) != 0:
#                 for record in tqdm(compounds_list, total=len(compounds_list), desc=f'Retrieving ligands for {accession}'):
#                     if record['molecule_structures']['molfile'] is not None:
#
#                         mol_id = record['molecule_chembl_id']
#                         file_name = f"{mol_id}_{accession}.sdf"
#                         file_path = os.path.join(save_path, "ligands", file_name)
#                         if os.path.exists(file_path):
#                             continue
#
#                         molfile = record['molecule_structures']['molfile']
#                         mw_freebase = f"> <mw_freebase>\n{record['molecule_properties']['mw_freebase']}\n"
#                         alogp = f"\n> <alogp>\n{record['molecule_properties']['alogp']}\n"
#                         mol_file += molfile
#                         mol_file += mw_freebase
#                         mol_file += alogp
#                         mol_file += "\n> <bioactivities>\n"
#
#                         pChEMBL_value = 0
#                         for i in bioactivities_list:
#                             if i['molecule_chembl_id'] == record['molecule_chembl_id']:
#                                 bioactivities_str = f"assay_type: {i['assay_type']}, pchembl_value: {i['pchembl_value']}, {i['standard_type']}{i['standard_relation']}{i['standard_value']}{i['standard_units']}\n"
#                                 mol_file += bioactivities_str
#                                 if i['pchembl_value'] is not None:
#                                     if float(i['pchembl_value']) > float(pChEMBL_value):
#                                         pChEMBL_value = str(i['pchembl_value'])
#                         if pChEMBL_value == 0 or pChEMBL_value == '0':
#                             pChEMBL_value = 'None'
#                         pChEMBL_value_str = f"\n> <pChEMBL_value>\n{pChEMBL_value}\n"
#                         mol_file += pChEMBL_value_str
#
#                         with open(file_path, 'w') as outfile:
#                             outfile.write(mol_file)
#                         mol_file = ''

from itertools import islice

def chunked_iterable(iterable, chunk_size):
    iterator = iter(iterable)
    while True:
        chunk = list(islice(iterator, chunk_size))
        if not chunk:
            break
        yield chunk

def download_chembl_compounds(accession, save_path, gene_name):
    # Create the "ligands" directory if it does not exist
    if not os.path.exists(os.path.join(save_path, "ligands")):
        os.makedirs(os.path.join(save_path, "ligands"))

    targets_api = new_client.target
    compounds_api = new_client.molecule
    bioactivities_api = new_client.activity

    targets = targets_api.get(
        target_components__accession=accession,
        target_type__in=['SINGLE PROTEIN', 'PROTEIN FAMILY']
    ).only("target_chembl_id")
    targets_id = list(set([target['target_chembl_id'] for target in targets]))

    if len(targets_id) != 0:
        target_columns_list = ['molecule_chembl_id', 'molecule_name', 'molecule_max_phase', 'molecular_weight', 'alogp',
                               'standard_type', 'standard_relation', 'standard_value', 'standard_units',
                               'pchembl_value', 'assay_type', 'target_name', 'target_organism', 'target_type']

        for target in targets_id:

            bioactivities = bioactivities_api.filter(target_chembl_id=target,
                                                     standard_type__in=["IC50", "EC50", "Ki", "Kd"],
                                                     assay_type__in=["B", "F"], relation__in=['=', '<', '<=']).only(
                *target_columns_list)

            bioactivities_list = list(bioactivities)

            compound_columns_list = ["molecule_chembl_id", 'molecule_structures', 'molecule_properties', 'mw_freebase']
            for bioactivity_chunk in chunked_iterable(bioactivities_list, 100):
                compounds = compounds_api.filter(
                    molecule_chembl_id__in=[x['molecule_chembl_id'] for x in bioactivity_chunk],
                    molecule_properties__mw_freebase__lte=700,
                    molfile__isnull=False
                ).only(*compound_columns_list)

                compounds_list = [record for record in compounds]

                if len(compounds_list) != 0:
                    for record in tqdm(compounds_list, total=len(compounds_list),
                                       desc=f'Retrieving ligands for {accession}'):
                        if record['molecule_structures']['molfile'] is not None:

                            mol_id = record['molecule_chembl_id']
                            file_name = f"{mol_id}_{accession}.sdf"
                            file_path = os.path.join(save_path, "ligands", file_name)
                            if os.path.exists(file_path):
                                continue
                            mol_file = ''
                            molfile = record['molecule_structures']['molfile']
                            mw_freebase = f"> <mw_freebase>\n{record['molecule_properties']['mw_freebase']}\n"
                            alogp = f"\n> <alogp>\n{record['molecule_properties']['alogp']}\n"
                            mol_file += molfile
                            mol_file += mw_freebase
                            mol_file += alogp
                            mol_file += "\n> <bioactivities>\n"

                            pChEMBL_value = 0
                            for i in bioactivities_list:
                                if i['molecule_chembl_id'] == record['molecule_chembl_id']:
                                    bioactivities_str = f"assay_type: {i['assay_type']}, pchembl_value: {i['pchembl_value']}, {i['standard_type']}{i['standard_relation']}{i['standard_value']}{i['standard_units']}\n"
                                    mol_file += bioactivities_str
                                    if i['pchembl_value'] is not None:
                                        if float(i['pchembl_value']) > float(pChEMBL_value):
                                            pChEMBL_value = str(i['pchembl_value'])
                            if pChEMBL_value == 0 or pChEMBL_value == '0':
                                pChEMBL_value = 'None'
                            pChEMBL_value_str = f"\n> <pChEMBL_value>\n{pChEMBL_value}\n"
                            mol_file += pChEMBL_value_str

                            with open(file_path, 'w') as outfile:
                                outfile.write(mol_file)