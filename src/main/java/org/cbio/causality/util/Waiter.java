package org.cbio.causality.util;

/**
 * @author Ozgun Babur
 */
public class Waiter
{
	public synchronized static void pause(long time)
	{
		Waiter w = new Waiter();
		w.bekle(time);
	}

	private synchronized void bekle(long time)
	{
		try
		{
			this.wait(time);
		}
		catch (InterruptedException e)
		{
			e.printStackTrace();
		}

	}
}
