package org.cbio.causality.binintanalysis;

import org.cbio.causality.data.portal.PortalDataset;
import org.cbio.causality.data.portal.PortalDatasetEnum;

/**
 * @author Ozgun Babur
 */
public enum Dataset1
{
	GBM(PortalDatasetEnum.GLIOBLASTOMA_MUT_CNCALL_EXPZ, PortalDatasetEnum.GLIOBLASTOMA_EXP),
	UCEC(PortalDatasetEnum.ENDOMETRIAL_MUT_CNCALL_EXPZ, PortalDatasetEnum.ENDOMETRIAL_EXP),
	BRCA(PortalDatasetEnum.BREAST_MUT_CNCALL_EXPZ, PortalDatasetEnum.BREAST_EXP),
	COADREAD(PortalDatasetEnum.COLON_MUT_CNCALL_EXPZ, PortalDatasetEnum.COLON_EXP),
	LUAD(PortalDatasetEnum.LUNG_MUT_CNCALL_EXPZ, PortalDatasetEnum.LUNG_EXP),
	THCA(PortalDatasetEnum.THYROID_MUT_CNCALL_EXPZ, PortalDatasetEnum.THYROID_EXP),
	LAML(PortalDatasetEnum.LEUKEMIA_MUT_CNCALL_EXPZ, PortalDatasetEnum.LEUKEMIA_EXP),
	SKCM(PortalDatasetEnum.MELANOMA_MUT_CNCALL_EXPZ, PortalDatasetEnum.MELANOMA_EXP);

	public PortalDataset mutCnCallExpZ;
	public PortalDataset exp;

	private Dataset1(PortalDatasetEnum mutCnCallExpZ, PortalDatasetEnum exp)
	{
		this.mutCnCallExpZ = mutCnCallExpZ;
		this.exp = exp;
	}

	public String code()
	{
		return name();
	}
}
