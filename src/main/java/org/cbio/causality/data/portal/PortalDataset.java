package org.cbio.causality.data.portal;


import java.lang.String;

/**
 * @author Ozgun Babur
 */
public enum PortalDataset
{
	GLIOBLASTOMA("gbm_tcga_pub", "gbm_tcga_pub_cnaseq", "gbm_tcga_pub_cna_consensus", "gbm_tcga_pub_mutations"),
	ENDOMETRIAL("ucec_tcga_pub", "ucec_tcga_pub_cnaseq", "ucec_tcga_pub_gistic", "ucec_tcga_pub_mutations");

	String cancerStudyID;
	String caseListID;
	String[] profileID;

	private PortalDataset(String cancerStudyID, String caseListID, String... profileID)
	{
		this.cancerStudyID = cancerStudyID;
		this.caseListID = caseListID;
		this.profileID = profileID;
	}
}
