package org.cbio.causality.data;

import org.cbio.causality.model.Alteration;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.model.Change;
import org.junit.Ignore;
import org.junit.Test;

/**
 * @author Ozgun Babur
 */
@Ignore //TODO it failed sometimes (NPE, file is null), sometimes - not... Should not depend on non-standard file path or internet connection (classpath locations can be used...)
public class GEOAccessorTest
{
	@Test
	public void readTheSeries()
	{
		String gseID = GSE11223.GSE_ID;
		String id = "657";
		int[] control = GSE11223.Normal_Uninflamed_sigmoid_colon;
		int[] test = GSE11223.UC_Uninflamed_sigmoid_colon;
		
		GEOAccessor acc = new GEOAccessor(gseID, test, control);
		AlterationPack pack = acc.getAlterations(id);
		for (Change ch : pack.get(Alteration.EXPRESSION))
		{
			System.out.println(ch);
		}
	}
}
