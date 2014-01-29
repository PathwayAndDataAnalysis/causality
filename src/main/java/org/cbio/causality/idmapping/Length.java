package org.cbio.causality.idmapping;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Gets length of the given gene.
 *
 * R code to get these length:
 *
 hg19GeneLengths <- function(symbols)
 {
	 require(TxDb.Hsapiens.UCSC.hg19.knownGene)
	 require(org.Hs.eg.db)
	 exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
	 egs = unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
	 sapply(egs,function(eg)
	 {
		 exons = exons.db[[eg]]
		 if (!is.null(exons))
		 {
			 exons = reduce(exons)
			 sum( width(exons) )
		 }
	 })
 }

 eg = read.table('/home/ozgun/Desktop/temp.txt')
 genes = as.vector(eg[,])
 result = hg19GeneLengths(genes)
 result = unlist(result)
 write.table(result, file="/home/ozgun/Desktop/gene-lengths.txt", sep="\t")

 *
 * @author Ozgun Babur
 */
public class Length
{
	private static Map<String, Integer> map;

	public static void main(String[] args)
	{
		System.out.println(Length.of("BAX"));
	}

	/**
	 * Provides HGNC ID of the given approved gene symbol.
	 * @param symbol
	 * @return
	 */
	public static Integer of(String symbol)
	{
		return map.get(symbol);
	}

	public static boolean containsSymbol(String symbol)
	{
		return map.containsKey(symbol);
	}

	public static Set<String> getSymbols()
	{
		return map.keySet();
	}

	public static int getTotalLength(Set<String> symbols)
	{
		int t = 0;
		for (String symbol : symbols)
		{
			Integer l = of(symbol);

			if (l != null) t += l;
		}
		return t;
	}

	static
	{
		try
		{
			map = new HashMap<String, Integer>();
			BufferedReader reader = new BufferedReader(new InputStreamReader(
				HGNC.class.getResourceAsStream("gene-lengths.txt")));
			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				String[] token = line.split("\t");
				String sym = token[0];
				Integer length = Integer.parseInt(token[1]);
				map.put(sym, length);
			}
			reader.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

}
