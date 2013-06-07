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
	OVARIAN_MUT("ov_tcga_pub", "brca_tcga_pub_sequenced", "brca_tcga_pub_mutations");

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
