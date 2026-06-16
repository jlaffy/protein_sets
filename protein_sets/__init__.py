"""Protein-set resolution: UniProt proteome queries, file-based sets, region tables, OMA orthologs."""

from .layout import (
    universe_path, attach_provenance,
    SUBGROUP_FROM_UNIPROT_FEATURE, SUBGROUP_FROM_INTERPRO_DB,
)
from .uniprot import (
    pull_proteome, load_proteome, PROTEOME_CACHE,
    pull_alt_isoforms, load_alt_isoforms,
)
from .uniref import (
    pull_uniref_mapping, load_uniref_mapping, fetch_cluster_ids, UNIREF_CACHE,
)
from .hpa import (
    pull_hpa_secretome, load_hpa_secretome, hpa_secretome_uniprots,
    pull_hpa_subloc, load_hpa_subloc, hpa_subloc,
    HPA_CACHE,
)
from .cspa import (
    pull_cspa, load_cspa_human, load_cspa_mouse, load_cspa_glycopeptides,
    cspa_uniprots, CSPA_CACHE,
)
from .humantfs import (
    pull_humantfs, load_humantfs, humantfs_uniprots, HUMANTFS_CACHE,
)
from .interpro import load_interpro, load_entry_types, INTERPRO_CACHE
from .oma import (
    OMA_CACHE,
    pull_uniprot_mapping, load_uniprot_mapping, best_accession_per_oma,
    pull_msa, pull_msas, load_msa, msa_path,
    load_coverage, pull_group_metadata, load_group_metadata,
    pull_hog_sizes, load_hog_sizes,
    oma_group_id, oma_msa, oma_msa_file, orthologs,
)
from .definitions import (
    reviewed_human, reviewed_human_alt_isoforms, all_isoforms_of,
    alt_isoforms_no_tm, alt_isoforms_lost_tm,
    signal_peptide, transmem,
    secretome, secretome_hpa, secretome_uniprot, secretome_with_sp,
    sp_strict_no_tm, no_sp_no_tm, no_secretome_no_tm,
    surfaceome, transcription_factors,
    NAMED_SETS, resolve_dynamic, list_dynamic,
    dedupe,
)
from . import studies
from .regions import (
    interpro_domain_regions, pfam_regions,
    from_external_table,
    uniprot_regions, uniprot_regions_remove,
    # backcompat aliases (importable, not advertised)
    domain_regions,
    uniprot_feature_regions, uniprot_feature_regions_remove,
    signal_peptide_regions, propeptide_regions,
    transit_peptide_regions, mito_transit_regions, mts_regions, pts2_regions,
    motif_regions, nls_regions, nes_regions, pts1_regions,
    er_retention_kdel_regions, er_retention_dibasic_regions,
    er_retention_regions, peroxisome_targeting_regions,
    mitochondrial_targeting_regions,
    nuclear_localization_regions, nuclear_localisation_regions,
    nuclear_export_regions,
    cell_attachment_rgd_regions, cell_attachment_atypical_regions,
    cell_attachment_regions, rgd_regions,
    pdz_binding_regions, sh3_binding_regions,
    lipidation_sites, gpi_sites, palmitoyl_sites, myristoyl_sites, prenyl_sites,
    glycosylation_sites, n_glycan_sites, o_glycan_sites, c_glycan_sites,
    disulfide_sites, disulfide_bonds,
    modified_residue_sites, phospho_sites, acetyl_sites, methyl_sites,
    UNIPROT_FEATURE_KEYWORDS,
)
from .lookup import resolve
from .subloc import uniprot_subloc, uniprot_subloc_wide
from .transforms import mutation_scan, to_fasta

__all__ = [
    'pull_proteome', 'load_proteome', 'PROTEOME_CACHE',
    'pull_alt_isoforms', 'load_alt_isoforms',
    'reviewed_human_alt_isoforms', 'all_isoforms_of',
    'alt_isoforms_no_tm', 'alt_isoforms_lost_tm',
    'pull_uniref_mapping', 'load_uniref_mapping', 'fetch_cluster_ids', 'UNIREF_CACHE',
    'pull_hpa_secretome', 'load_hpa_secretome', 'hpa_secretome_uniprots',
    'pull_hpa_subloc', 'load_hpa_subloc', 'hpa_subloc',
    'HPA_CACHE',
    'pull_cspa', 'load_cspa_human', 'load_cspa_mouse', 'load_cspa_glycopeptides',
    'cspa_uniprots', 'CSPA_CACHE',
    'pull_humantfs', 'load_humantfs', 'humantfs_uniprots', 'HUMANTFS_CACHE',
    'load_interpro', 'load_entry_types', 'INTERPRO_CACHE',
    'OMA_CACHE',
    'pull_uniprot_mapping', 'load_uniprot_mapping', 'best_accession_per_oma',
    'pull_msa', 'pull_msas', 'load_msa', 'msa_path',
    'load_coverage', 'pull_group_metadata', 'load_group_metadata',
    'pull_hog_sizes', 'load_hog_sizes',
    'oma_group_id', 'oma_msa', 'oma_msa_file', 'orthologs',
    'reviewed_human', 'signal_peptide', 'transmem',
    'secretome', 'secretome_hpa', 'secretome_uniprot', 'secretome_with_sp',
    'sp_strict_no_tm', 'no_sp_no_tm', 'no_secretome_no_tm',
    'surfaceome', 'transcription_factors',
    'NAMED_SETS', 'resolve_dynamic', 'list_dynamic',
    'dedupe',
    'studies',
    'interpro_domain_regions', 'pfam_regions',
    'from_external_table',
    'uniprot_regions', 'uniprot_regions_remove',
    'domain_regions',
    'uniprot_feature_regions', 'uniprot_feature_regions_remove',
    'signal_peptide_regions', 'propeptide_regions',
    'transit_peptide_regions', 'mito_transit_regions', 'mts_regions', 'pts2_regions',
    'motif_regions', 'nls_regions', 'nes_regions', 'pts1_regions',
    'er_retention_kdel_regions', 'er_retention_dibasic_regions',
    'er_retention_regions', 'peroxisome_targeting_regions',
    'mitochondrial_targeting_regions',
    'nuclear_localization_regions', 'nuclear_localisation_regions',
    'nuclear_export_regions',
    'cell_attachment_rgd_regions', 'cell_attachment_atypical_regions',
    'cell_attachment_regions', 'rgd_regions',
    'pdz_binding_regions', 'sh3_binding_regions',
    'lipidation_sites', 'gpi_sites', 'palmitoyl_sites', 'myristoyl_sites', 'prenyl_sites',
    'glycosylation_sites', 'n_glycan_sites', 'o_glycan_sites', 'c_glycan_sites',
    'disulfide_sites', 'disulfide_bonds',
    'modified_residue_sites', 'phospho_sites', 'acetyl_sites', 'methyl_sites',
    'UNIPROT_FEATURE_KEYWORDS',
    'resolve',
    'uniprot_subloc', 'uniprot_subloc_wide',
    'mutation_scan', 'to_fasta',
]
