package org.cbio.causality.data.portal;


import java.lang.String;

/**
 * @author Ozgun Babur
 */
public enum PortalDataset
{
	GLIOBLASTOMA_MUT_CN("gbm_tcga_pub", "gbm_tcga_pub_cnaseq", "gbm_tcga_pub_cna_consensus", "gbm_tcga_pub_mutations"),
	ENDOMETRIAL_MUT_CN("ucec_tcga_pub", "ucec_tcga_pub_cnaseq", "ucec_tcga_pub_gistic", "ucec_tcga_pub_mutations"),
	ENDOMETRIAL_MUT("ucec_tcga_pub", "ucec_tcga_pub_sequenced", "ucec_tcga_pub_mutations"),
	ENDOMETRIAL_MUT_HYPER("ucec_tcga_pub", "ucec_tcga_pub_msi", "ucec_tcga_pub_mutations"),
	ENDOMETRIAL_MUT_ULTRA("ucec_tcga_pub", "ucec_tcga_pub_pole", "ucec_tcga_pub_mutations"),
	BREAST_MUT("brca_tcga_pub", "brca_tcga_pub_sequenced", "brca_tcga_pub_mutations"),
	BREAST_MUT_CN("brca_tcga_pub", "brca_tcga_pub_cnaseq", "brca_tcga_pub_mutations", "brca_tcga_pub_gistic"),
	LUNG_ADEN_MUT_CN("luad_tcga", "luad_tcga_cnaseq", "luad_tcga_mutations", "luad_tcga_gistic"),
	COLON_MUT_CN("coadread_tcga_pub", "coadread_tcga_pub_cna_seq", "coadread_tcga_pub_mutations", "coadread_tcga_pub_gistic"),
	PROSTATE_MSKCC_MUT_CN("prad_mskcc", "prad_mskcc_cna_seq", "prad_mskcc_mutations", "prad_mskcc_cna"),
	PROSTATE_TCGA_MUT_CN("prad_tcga", "prad_tcga_cnaseq", "prad_tcga_mutations", "prad_tcga_gistic"),
	KIDNEY_MUT_CN("kirc_tcga_pub", "kirc_tcga_pub_cnaseq", "kirc_tcga_pub_mutations", "kirc_tcga_pub_gistic"),
	OVARIAN_MUT("ov_tcga_pub", "ov_tcga_pub_sequenced", "ov_tcga_pub_mutations"),
	OVARIAN_MUT_CN("ov_tcga_pub", "ov_tcga_pub_cna_seq", "ov_tcga_pub_mutations", "ov_tcga_pub_gistic");

	public String cancerStudyID;
	public String caseListID;
	public String[] profileID;

	private PortalDataset(String cancerStudyID, String caseListID, String... profileID)
	{
		this.cancerStudyID = cancerStudyID;
		this.caseListID = caseListID;
		this.profileID = profileID;
	}
}
