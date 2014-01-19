package org.cbio.causality.util;

/**
 * @author Ozgun Babur
 */
public class FormatUtil
{
	public static double roundToSignificantDigits(double v, int digits)
	{
		double x = v;

		int a = 0;
		while (x < 1)
		{
			x *= 10;
			a++;
		}

		if (a == 0)
		{
			while (x > 1)
			{
				x /= 10;
				a--;
			}
		}

		int shift = a + digits;

		double c = 1;

		if (shift > 0)
		{
			for (int i = 0; i < shift; i++) c *= 10;
		}
		else if (shift < 0)
		{
			for (int i = 0; i < -shift; i++) c /= 10;
		}

		return Math.round(v * c) / c;
	}
}
