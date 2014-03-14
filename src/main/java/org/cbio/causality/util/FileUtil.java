package org.cbio.causality.util;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;

import java.io.*;
import java.util.Enumeration;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * @author Ozgun Babur
 */
public class FileUtil
{
	public static ZipEntry findEntryContainingNameInZIPFile(String zipFileName,
		String partOfEntryName)
	{
		try
		{
			ZipFile zipFile = new ZipFile(zipFileName);
			Enumeration<? extends ZipEntry> entries = zipFile.entries();

			while (entries.hasMoreElements())
			{
				ZipEntry zipEntry = entries.nextElement();
				if (!zipEntry.isDirectory())
				{
					String fileName = zipEntry.getName();
					if (fileName.contains(partOfEntryName))
					{
						return zipEntry;
					}
				}
			}
			zipFile.close();
		} catch (final IOException ioe)
		{
			ioe.printStackTrace();
			return null;
		}
		return null;
	}

	public static boolean extractEntryContainingNameInTARGZFile(String targzFileName,
		String partOfEntryName, String extractedName)
	{
		try
		{
			TarArchiveEntry entry;

			TarArchiveInputStream is = new TarArchiveInputStream(new GZIPInputStream(
				new FileInputStream(targzFileName)));

			while ((entry = is.getNextTarEntry()) != null)
			{
				if (entry.isDirectory()) continue;
				else
				{
					if (entry.getName().contains(partOfEntryName))
					{
						byte [] btoRead = new byte[1024];

						BufferedOutputStream bout =new BufferedOutputStream(
							new FileOutputStream(extractedName));

						int len;
						while((len = is.read(btoRead)) != -1)
						{
							bout.write(btoRead,0,len);
						}

						bout.close();

						return true;

					}
				}
			}
		}
		catch (IOException ioe)
		{
			ioe.printStackTrace();
		}

		return false;
	}

	public static String getFileContent(String filename)
	{
		try
		{
			StringBuilder sb = new StringBuilder();
			BufferedReader reader = new BufferedReader(new FileReader(filename));
			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				sb.append(line).append("\n");
			}

			reader.close();
			return sb.toString();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return null;
	}

	public static void printLines(String filename, String partialContent) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));

		int i = 0;
		for (String line = reader.readLine(); line != null; line = reader.readLine())
		{
			i++;
			if (line.contains(partialContent))
			{
				System.out.println("Line " + i + ": " + line);
			}
		}

		reader.close();
	}

	public static void main(String[] args) throws IOException
	{
		printLines("SIFWithLoc.sif", "controls-transport-of");
	}
}
