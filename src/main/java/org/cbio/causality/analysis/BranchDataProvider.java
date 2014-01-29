package org.cbio.causality.analysis;

import java.awt.*;

/**
 * @author Ozgun Babur
 */
public interface BranchDataProvider
{
	public Color getColor(String gene, String root);
	public double getThickness(GeneBranch branch, String root);
}
