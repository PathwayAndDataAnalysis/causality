package org.cbio.causality.binintanalysis;

import org.cbio.causality.data.portal.PortalDataset;

/**
 * @author Ozgun Babur
 */
public enum Dataset1
{
	GBM(PortalDataset.GLIOBLASTOMA_MUT_CNCALL_EXPZ, PortalDataset.GLIOBLASTOMA_EXP),
	UCEC(PortalDataset.ENDOMETRIAL_MUT_CNCALL_EXPZ, PortalDataset.ENDOMETRIAL_EXP),
	BRCA(PortalDataset.BREAST_MUT_CNCALL_EXPZ, PortalDataset.BREAST_EXP),
	COAD(PortalDataset.COLON_MUT_CNCALL_EXPZ, PortalDataset.COLON_EXP),
	LUAD(PortalDataset.LUNG_MUT_CNCALL_EXPZ, PortalDataset.LUNG_EXP),
	THCA(PortalDataset.THYROID_MUT_CNCALL_EXPZ, PortalDataset.THYROID_EXP),
	LAML(PortalDataset.LEUKEMIA_MUT_CNCALL_EXPZ, PortalDataset.LEUKEMIA_EXP);

	PortalDataset mutCnCallExpZ;
	PortalDataset exp;

	private Dataset1(PortalDataset mutCnCallExpZ, PortalDataset exp)
	{
		this.mutCnCallExpZ = mutCnCallExpZ;
		this.exp = exp;
	}

	public String code()
	{
		return name();
	}
}
