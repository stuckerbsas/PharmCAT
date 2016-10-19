package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import javax.annotation.Nonnull;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.w3c.dom.Document;


/**
 * This extracts the positions from the definition files and reformats them as a vcf.
 * As a command line it takes three arguments:
 * -d  CPIC definitions directory
 * -o output vcf
 * -g genome build to use for DAS
 *
 * Created by lester on 6/20/16.
 */

public class ExtractPositions {
  private static Path ps_definitionDir;
  private static Path ps_outputVcf;
  private static String m_genomeBuild;
  private static final String sf_fileHeader = "##fileformat=VCFv4.1\n" +
      "##fileDate=2015-08-04\n" +
      "##source=Electronic, version: hg38_2.0.1\n" +
      "##reference=hg38\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "#CHROM  POS  ID REF  ALT  QUAL FILTER INFO FORMAT PharmCAT\n";


  // Default constructor
  public ExtractPositions(Path definitionDir, Path outputVcf, String genomeBuild) {
    ps_definitionDir = definitionDir;
    ps_outputVcf = outputVcf;
    m_genomeBuild = genomeBuild;
  }


  // Main method for command line use
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("d", "definition-dir", "directory of allele definition files", true, "d")
          .addOption("o", "output-file", "output vcf file", true, "o")
          .addOption("g", "genome-build", "hg38 or hg37", true, "g");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path definitionDir = cliHelper.getValidDirectory("d", false);
      Path outputVcf= cliHelper.getValidFile("o", false);

      if (!cliHelper.getValue("g").equalsIgnoreCase("hg38") && !cliHelper.getValue("g").equalsIgnoreCase("hg37")) {
        System.out.println("-g parameter must be hg37 or hg38");
        System.exit(1);
      }
      String genomeBuild = cliHelper.getValue("g");

      ExtractPositions extractPositions = new ExtractPositions(definitionDir, outputVcf, genomeBuild);
      extractPositions.run();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  //Extract the positions and print them to file
  private void run() {
    try {
      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(ps_definitionDir);
      if (definitionReader.getGenes().size() == 0) {
        System.out.println("Did not find any allele definitions at " + ps_definitionDir);
        System.exit(1);
      }
      StringBuilder vcfString = getPositions(definitionReader);
      try (PrintWriter writer = new PrintWriter(String.valueOf(ps_outputVcf), "UTF-8")) {  // close PrintWriter with try
        writer.print(vcfString);
        writer.flush();
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }

  }


  // Build up vcf string
  public  StringBuilder getPositions(@Nonnull DefinitionReader definitionReader) throws IOException {
    StringBuilder builder = new StringBuilder();
    builder.append(sf_fileHeader);
    for (String gene : definitionReader.getGenes()) {  // for each definition file
      for (VariantLocus variantLocus :definitionReader.getPositions(gene)) {
          String[] vcfFields = getVcfLineFromDefinition(definitionReader, gene, variantLocus, m_genomeBuild);
          String vcfLine = String.join("\t", (CharSequence[])vcfFields);
          builder.append(vcfLine+"\n");
      }
    }
    return  builder;
  }

  /*
  * Helper method to convert repeats into standard format
  */
  public static String expandAllele(@Nonnull String allele) {
    String finalAllele = allele;
    if (allele.contains("(")) {
      int bracketStart=0;
      int bracketEnd=0;
      String repeat = "";
      for (int i = 0; i < allele.toCharArray().length; i++) {
        String character = String.valueOf(allele.toCharArray()[i]);
        if (character.equals("(")) {
          bracketStart = i;
        }
        if (character.equals(")")) {
          bracketEnd = i;
        }
        if (character.matches("[0-9]+")) {
          repeat = repeat + character; // just in case we have a repeat greater than nine

        }
      }
      String repeatSeq = allele.substring(bracketStart+1, bracketEnd);
      String expandedAlelle= "";
      for (int repeats =0; repeats< Integer.parseInt(repeat); repeats++) {
        expandedAlelle = expandedAlelle + repeatSeq;
      }
      finalAllele = allele.substring(0, bracketStart) + expandedAlelle + allele.substring(bracketEnd+1+repeat.length(), allele.length());
    }
    return finalAllele;
  }


  /*
  * For deletions we need to be able to get the nucleotide at the previous position.
  * This is not stored in the definition file, so we need to use an external source.
  *
   */
  public static String getDAS(String chr, String position, String genomeBuild) {
    String nucleotide= "";
    try {
      String uri =
          "http://genome.ucsc.edu/cgi-bin/das/"+ genomeBuild +  "/dna?segment="+ chr +":" +position +"," + position;
      System.out.println("Getting position: " + uri);
      URL url = new URL(uri);
      HttpURLConnection connection =
          (HttpURLConnection) url.openConnection();
      connection.setRequestMethod("GET");
      connection.setRequestProperty("Accept", "application/xml");

      InputStream xml = connection.getInputStream();

      DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
      DocumentBuilder db = dbf.newDocumentBuilder();
      Document doc = db.parse(xml);
      doc.getDocumentElement().normalize();

      nucleotide = doc.getElementsByTagName("DNA").item(0).getTextContent();
    }  catch (Exception e) {
      e.printStackTrace();
    }
    nucleotide = nucleotide.toUpperCase().trim();
  return nucleotide;
  }


  /*
  * Helper method to convert repeats into standard format
  */
  public static String[] getVcfLineFromDefinition(@Nonnull DefinitionReader definitionReader, @Nonnull String gene, @Nonnull VariantLocus variantLocus,
      String genomeBuild) {
    DefinitionFile definitionFile = definitionReader.getDefinitionFile(gene);
    NamedAllele namedAlleleFirst = definitionFile.getNamedAlleles().get(0); //Get star one for Ref column of vcf TODO:: use DAS?
    String allele = namedAlleleFirst.getAllele(variantLocus);
    int position = variantLocus.getVcfPosition();
    String rsid = variantLocus.getRsid();
    if (rsid == null) {
      rsid = ".";
    }

    String chr = definitionReader.getDefinitionFile(gene).getChromosome();
    ArrayList<String> alts = new ArrayList<>(); //get alts from the namedAlleles
    ArrayList<String> starAlleles = new ArrayList<>();
    definitionFile.getNamedAlleles().stream().filter(namedAllele -> namedAllele.getAllele(variantLocus) != null
      ).forEachOrdered(namedAllele -> {
      if (!alts.contains(expandAllele(namedAllele.getAllele(variantLocus)))
          && !allele.equals(expandAllele(namedAllele.getAllele(variantLocus)))) {
        alts.add(expandAllele(namedAllele.getAllele(variantLocus)));
      }

      //calculate size of haplotype
      ArrayList l = new ArrayList(Arrays.asList(namedAllele.getAlleles()));
      l.removeAll(Collections.singleton(null));


      starAlleles.add(gene + ":" + namedAllele.getName().replace(" ","")+"[" + l.size() +"]" +"=" + namedAllele.getAllele(variantLocus));
    });

    if (alts.size() == 0 ) {
      alts.add(".");
    }

    String finalAllele = expandAllele(allele);

    if (alts.contains("Y")) {
      alts.remove("Y");
      if (!alts.contains("C") && !allele.equals("C")) {
        alts.add("C");
      }
      if (!alts.contains("T")  && !allele.equals("T")) {
        alts.add("T");
      }
    }
    if (alts.contains("R")) {
      alts.remove("R");
      if (!alts.contains("G")  && !allele.equals("G")) {
        alts.add("G");
      }
      if (!alts.contains("A")  && !allele.equals("A")) {
        alts.add("A");
      }
    }

    // Simplest possible del/ins parsing, presuming the only a single ins or del string, or the correct one being first
    if (alts.size()>0) {
      if (alts.get(0).contains("ins")) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        alts.set(0, nucleotide + alts.get(0).replace("ins", ""));
        finalAllele = nucleotide;
      }
      if (alts.get(0).contains("del")) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        finalAllele = nucleotide+ alts.get(0).replace("del", "");
        alts.set(0, nucleotide);
      }
    }

    String alt = String.join(",", alts);
    String starAllele = String.join(",", starAlleles);

    String[] vcfFields = {
        chr,
        Integer.toString(position),
        rsid,
        finalAllele,
        alt,
        ".",
        "PASS",
        starAllele,
        "GT",
        "0/0" };
    return vcfFields ;
    }
}

