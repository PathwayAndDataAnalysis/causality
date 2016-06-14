package org.cbio.causality.util;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;

import java.io.*;
import java.nio.channels.FileChannel;
import java.util.Enumeration;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;

/**
 * @author Ozgun Babur
 */
public class FileUtil
{
	public static ZipEntry findEntryContainingNameInZIPFile(String zipFileName,
		String partOfEntryName)
	{
		return findEntryContainingNameInZIPFile(zipFileName, partOfEntryName, null);
	}

	public static ZipEntry findEntryContainingNameInZIPFile(String zipFileName,
		String partOfEntryName, String tabooPartOfEntryName)
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
					if (fileName.contains(partOfEntryName) &&
						(tabooPartOfEntryName == null || !fileName.contains(tabooPartOfEntryName)))
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
						is.close();

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

	public static boolean extractAllEntriesContainingNameInTARGZFile(String targzFileName,
		String partOfEntryName, String dir)
	{
		try
		{
			TarArchiveEntry entry;

			TarArchiveInputStream is = new TarArchiveInputStream(new GZIPInputStream(
				new FileInputStream(targzFileName)));

			boolean success = false;

			while ((entry = is.getNextTarEntry()) != null)
			{
				if (entry.isDirectory()) continue;
				else
				{
					String name = entry.getName();
					if (name.contains("/")) name = name.substring(name.lastIndexOf("/") + 1);
					if (name.contains(partOfEntryName))
					{
						byte [] btoRead = new byte[1024];

						BufferedOutputStream bout =new BufferedOutputStream(
							new FileOutputStream(dir + File.separator + name));

						int len;
						while((len = is.read(btoRead)) != -1)
						{
							bout.write(btoRead,0,len);
						}

						bout.close();
						success = true;
					}
				}
			}
			is.close();
			return success;
		}
		catch (IOException ioe)
		{
			ioe.printStackTrace();
		}

		return false;
	}

	public static boolean extractEntryContainingNameInZipFile(String zipFileName,
		String partOfEntryName, String tabooPartOfEntryName, String extractedName)
	{
		try
		{
			ZipEntry entry = findEntryContainingNameInZIPFile(zipFileName, partOfEntryName,
				tabooPartOfEntryName);

			OutputStream out = new FileOutputStream(extractedName);
			FileInputStream fin = new FileInputStream(zipFileName);
			BufferedInputStream bin = new BufferedInputStream(fin);
			ZipInputStream zin = new ZipInputStream(bin);
			ZipEntry ze;
			while ((ze = zin.getNextEntry()) != null) {
				if (ze.getName().equals(entry.getName())) {
					byte[] buffer = new byte[8192];
					int len;
					while ((len = zin.read(buffer)) != -1) {
						out.write(buffer, 0, len);
					}
					out.close();
					break;
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

	public static void printLines(String filename, int fromLine, int toLine) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));

		int i = 0;
		for (String line = reader.readLine(); line != null; line = reader.readLine())
		{
			i++;
			if (i >= fromLine && i <= toLine)
			{
				System.out.println(line);
			}
			if (i > toLine) break;
		}

		reader.close();
	}

	public static void copyFile(String src, String dest) throws IOException
	{
		File sourceFile = new File(src);
		File destFile = new File(dest);

		if(!destFile.exists())
		{
			destFile.createNewFile();
		}

		FileChannel source = null;
		FileChannel destination = null;

		try {
			source = new FileInputStream(sourceFile).getChannel();
			destination = new FileOutputStream(destFile).getChannel();
			destination.transferFrom(source, 0, source.size());
		}
		catch (IOException e){e.printStackTrace();}
		finally {
			if(source != null) {
				source.close();
			}
			if(destination != null) {
				destination.close();
			}
		}
	}

	public static void delete(File dir)
	{
		if (dir.isDirectory())
		{
			for (File file : dir.listFiles())
			{
				delete(file);
			}
		}
		dir.delete();
	}

	public static void main(String[] args) throws IOException
	{
//		printLines("SIFWithLoc.sif", "controls-transport-of");
		printLines("/home/babur/Projects/causalpath/src/main/resources/org/babur/causalpath/resource/mutation-stats.txt", 1, 1);
	}
}
