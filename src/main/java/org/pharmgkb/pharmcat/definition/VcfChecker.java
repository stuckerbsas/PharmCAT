package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;


/**
 *
 * Checks vcf for format issues in the reference and alt alleles
 *
 * Created by lester on 6/23/16.
 *
 * The following command line options are needed:
 * -d directory of allele definition files
 * -i input vcf file
 * -o output cleaned vcf file
 * -l look for liftover positions
 * -p run PharmCAT as well
 * -m add missing positions from the definition files
 * -g genome build, default is grc38
 */

public class VcfChecker {
  private  Path m_inputVcf;
  private  Path m_outputVcf;
  private boolean m_pharmcat;
  private boolean m_addMissing;
  private static String m_genomeBuild;
  private  LinkedHashMap<String, ArrayList<String>> m_vcfMap = new LinkedHashMap<>();
  private  HashMap<String, ArrayList<String>> m_vcfDefinitionLine = new HashMap<>();
  private  LinkedHashMap<String, ArrayList<String>> m_vcfMapFixed = new LinkedHashMap<>();
  private LinkedHashMap<String, String> m_vcfMapRepairs = new LinkedHashMap<>();
  private  DefinitionReader m_definitionReader = new DefinitionReader();
  private  ArrayList<String> m_vcfHeader = new ArrayList<>();
  private StringBuilder m_pharmCATResults = new StringBuilder();
  private StringBuilder m_pharmCATSummary = new StringBuilder();


  // Default constructor
  public VcfChecker(Path definitionDir, Path inputVcf, Path outputVcf, boolean pharmcat, String genomeBuild, boolean addMissing) {
    m_inputVcf = inputVcf;
    m_outputVcf = outputVcf;
    m_genomeBuild = genomeBuild;
    m_pharmcat = pharmcat;
    m_addMissing = addMissing;
    try {
      m_definitionReader.read(definitionDir);
      if (m_definitionReader.getGenes().size() == 0) {
        System.out.println("Did not find any allele definitions at " + definitionDir);
        System.exit(1);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  // Main method for command line use
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("d", "definition-dir", "directory of allele definition files", true, "d")
          .addOption("i", "input-vcf", "input vcf file", true, "i")
          .addOption("o", "output-vcf", "output vcf file", true, "o")
          .addOption("p", "pharmcat", "Run PharmCAT on fixed file")
          .addOption("m", "add-missing", "Add missing positions as star one")
          .addOption("g", "genome-build", "genome build hg37 or default hg38", false, "g");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path definitionDir = cliHelper.getValidDirectory("d", false);
      Path inputVcf= cliHelper.getValidFile("i", false);
      Path outputVcf= cliHelper.getValidFile("o", false);

      boolean pharmcat = false;
      if (cliHelper.hasOption("p")) {
        pharmcat = true;
      }

      boolean addMissing = false;
      if (cliHelper.hasOption("m")) {
        addMissing = true;
      }

      String genomeBuild = "hg38";
      if (cliHelper.hasOption("g")) {
        if ( !cliHelper.getValue("g").equalsIgnoreCase("hg38")) {
          System.out.println("-g parameter must present and be hg37 or hg38");
          System.exit(1);
        }
        genomeBuild = cliHelper.getValue("g");
      }

      VcfChecker vcfChecker = new VcfChecker(definitionDir, inputVcf, outputVcf, pharmcat, genomeBuild, addMissing);
      vcfChecker.run();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  //Extract the positions and print them to file
  public void run() {
      parseVcf();
      addDefinitionVcfLines();
      checkVcf();
      checkMissingVcfLines();
      writeVcf();
      if (m_pharmcat) {
        runPharmCat();
      }
      vcfToHtml();

  }


  private void runPharmCat() {
    NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(m_definitionReader, true, false);

    m_pharmCATResults.append("<H3>PharmCAT Results:</H3>");
    m_pharmCATResults.append("<table id=\"example2\" class=\"table table-bordered\" cellspacing=\"0\" width=\"520px\">\n" +
        "<thead><tr> <th>Gene</th> <th>Calls</th> <th>Missing</th> <th>Uncallable</th></tr></thead>\n<tbody>");
    m_pharmCATSummary.append("Overall Summary:" );
    m_pharmCATSummary.append(m_outputVcf).append("\t");

    try {
      System.out.println("Calling file: " + m_outputVcf);
      Result result = namedAlleleMatcher.call(m_outputVcf);
      new ResultSerializer()
          .alwaysShowUnmatchedHaplotypes(true)
          .toHtml(result, Paths.get(m_outputVcf.toString()+".pc.html"));
      for (GeneCall geneCalls : result.getGeneCalls()) {
        System.out.println("\nGene: " + geneCalls.getGene());
        ArrayList<String> diplotypes = new ArrayList<>();
        geneCalls.getDiplotypes().forEach(f-> diplotypes.add (f.toString() +" (score:"+ f.getScore() + ") "));
        System.out.print("Diplotypes ("+ geneCalls.getGene() +"): ");
        System.out.println(String.join( " ", diplotypes));
        ArrayList<String> missing = new ArrayList<>();
        missing.add("Missing positions  (" + geneCalls.getMatchData().getMissingPositions().size() + " out of " +
            (geneCalls.getMatchData().getMissingPositions().size() +  geneCalls.getVariants().size()) + "): ");
        geneCalls.getMatchData().getMissingPositions().forEach(f-> missing.add(geneCalls.getChromosome() + " " + f.getVcfPosition() + "   "));
        System.out.println(String.join(" ", missing));

        ArrayList<String> uncallable = new ArrayList<>();
        uncallable.add("Uncallable (" + geneCalls.getUncallableHaplotypes().size()  + "): " );
        geneCalls.getUncallableHaplotypes().forEach(f->uncallable.add(f + "  "));
        System.out.println(String.join(" ", uncallable));
        m_pharmCATSummary.append(geneCalls.getGene()).append(":").append(String.join(" ", diplotypes)).append("\t");

        m_pharmCATResults.append("<tr> <td><a href=\"").append(m_outputVcf.toString()).append(".pc.html").append("#")
            .append(geneCalls.getGene()).append("\"> ").append(geneCalls.getGene()).append("</a></td> <td>")
            .append(String.join(" ", diplotypes)).append("</td> <td>").append(String.join(" ", missing)).append("</td> <td>")
            .append(String.join(" ", uncallable)).append("</td></tr> ");
      }
      m_pharmCATResults.append(" \n</tbody></table>");
      System.out.println(m_pharmCATSummary);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }


  private void checkVcf() {
    /*
    This interactively checks the vcf, variant by variant.
     */
    Scanner s = new Scanner(System.in);
    for (Map.Entry<String,ArrayList<String>> entry : m_vcfMap.entrySet()) {
      System.out.print("Checking line:" + String.join("\t", entry.getValue()));
      ArrayList<String> vcfLineFromDefinition = getVcfLineFromDefinition(entry.getKey());

      boolean found = checkPresent(entry.getKey());
      if (!found) {
        System.out.println("---> Line is not present in definition files ");
      }
      else {
        String chr = entry.getValue().get(0);
        String genotype = getGenotype(entry.getValue().get(8), entry.getValue().get(9));

        //update rsids to our rsids
        if (getVcfLineFromDefinition(entry.getKey()).get(2) != null) {
          entry.getValue().set(2, getVcfLineFromDefinition(entry.getKey()).get(2));
        }

        // 1) Check if vcf variants start with chr. If not suggest they fix using grep, as it might be difficult to fix.
        if (!chr.startsWith("chr")) {
          System.out.println("Chromosomes must be in chr1, chr22, etc format. Please fix the the file and start again");
          System.out.println("Your can try using this command to fix this:");
          System.out.println("");
          System.exit(1);
        }

        //2) Check liftover - naive flip of genotypes from NCBI REMAP
        if (entry.getValue().get(7).contains("REF_UPDATE")) {
          System.out.println("---> This position has been lifted over. Updating genotype.");
          String newGenotype = genotype.replace("0", "N");
          newGenotype = newGenotype.replace("1", "0");
          newGenotype = newGenotype.replace("N", "1");
          recordRepair(entry.getKey(), "Liftover - genotype updated from " + genotype + " to " + newGenotype);
          entry.getValue().set(9, entry.getValue().get(9).replace(genotype, newGenotype));
          entry.getValue().set(7, vcfLineFromDefinition.get(7)+";"+ "LIFTOVER_GENOTYPE_SWITCHED;" + entry.getValue().get(7));
        }

        // 3) Line is present, but as it is not readable (./.) it can be ignored
        if (genotype.equals("./.") || genotype.equals(".|.)")) {
          System.out.println("---> Line is present, but is ./. or .|. so will be ignored. Added.");
          recordRepair(entry.getKey(), "Line is present, but is ./. or .|.");
          entry.getValue().set(7, vcfLineFromDefinition.get(7)+";"+ "UNREADABLE;");
          m_vcfMapFixed.put(entry.getKey(), entry.getValue());
        } else {
          // 3) Ref and Alts match! - can add line
          assert vcfLineFromDefinition != null;
          if (entry.getValue().get(3).equals(vcfLineFromDefinition.get(3)) && entry.getValue().get(4).equals(vcfLineFromDefinition.get(4))) {
            System.out.println("---> Ref and Alt positions match definition file.  Added");
            recordRepair(entry.getKey(), "Line is present, ref and alt matched.");
            entry.getValue().set(7, vcfLineFromDefinition.get(7)+";"+ "CORRECT;" + entry.getValue().get(7));
            m_vcfMapFixed.put(entry.getKey(), entry.getValue());
          }
            // 4) Alt does not match, but as genotype is 0/0 it doesn't matter. Likely no call
            else if (entry.getValue().get(3).equals(vcfLineFromDefinition.get(3)) && !entry.getValue().get(4).equals(vcfLineFromDefinition.get(4))
              && (genotype.equals("0/0") || genotype.equals("0|0"))) {
            System.out.println("---> Alt does not match, but genotype is 0/0 or 0|0 so doesn't matter. Added.");
            recordRepair(entry.getKey(), "Alt is incorrect (should be:" + vcfLineFromDefinition.get(4)+"), but is 0/0 or 0|0 so doesn't matter");
            entry.getValue().set(7, vcfLineFromDefinition.get(7)+";"+ "WRONG_ALT_IRRELEVANT;" + entry.getValue().get(7));
            m_vcfMapFixed.put(entry.getKey(), entry.getValue());
          }
            // 5) Ref and alt don't match. Update to Definition file values, or change.
            else {
            System.out.println("---> Ref or Alt does not match expected (" +String.join("\t", vcfLineFromDefinition).trim() + "). Update line to:");
            recordRepair(entry.getKey(),"Repaired (ref was :" + entry.getValue().get(3) +
                " alt was:"+ entry.getValue().get(4) +" )");
            entry.getValue().set(7, vcfLineFromDefinition.get(7)+";"+ "UPDATED_FROM_" + entry.getValue().get(3) +"_"+ entry.getValue().get(4)+";" + entry.getValue().get(7));

            entry.getValue().set(3, vcfLineFromDefinition.get(3));
            entry.getValue().set(4, vcfLineFromDefinition.get(4));
            System.out.println("------> "+String.join("\t", entry.getValue()));
            System.out.println("------> Does this line work (default is yes):");
            String  updateDecision = s.nextLine();
            if (!updateDecision.isEmpty() && !updateDecision.equalsIgnoreCase("y") && !updateDecision.equalsIgnoreCase("yes")) {
              System.out.println("Type in new ref:");
              String newRef = s.nextLine();
              System.out.println("Type in new alt:");
              String newAlt = s.nextLine();
              entry.getValue().set(3, newRef);
              entry.getValue().set(4, newAlt);
              System.out.println("------> Adding new line: " + String.join("\t", entry.getValue()));
              m_vcfMapFixed.put(entry.getKey(), entry.getValue());
            } else {
              System.out.println("Line added");
              m_vcfMapFixed.put(entry.getKey(), entry.getValue());
            }
          }
        }
      }
    }
  }


  // adds the repair text to the repairs ArrayList
  private void recordRepair(String key, String repairNote) {
    if (m_vcfMapRepairs.containsKey(key)) {
      m_vcfMapRepairs.put(key, m_vcfMapRepairs.get(key) + " " + repairNote);
    } else {
      m_vcfMapRepairs.put(key, repairNote);

    }
  }


  // Get the vcf line
  private ArrayList<String> getVcfLineFromDefinition(String position) {
    if (m_vcfDefinitionLine.containsKey(position)) {
      return m_vcfDefinitionLine.get(position);
    } else {
      return null;
    }
  }


  //Get the genotype
  private static String getGenotype(String infoString, String dataString) {
    String info[] = infoString.split(":");
    String data[] = dataString.split(":");

    String genotype = "";
    for(int i=0; i< info.length; i++){
      if (info[i].equals("GT") ){
        genotype = data[i];
      }
    }
    return genotype;
  }


  //Gets scores like DP from format and records
  private static String getScores(String infoString, String dataString) {
    String info[] = infoString.split(":");
    String data[] = dataString.split(":");

    String scores= "";
    for(int i=0; i< info.length; i++){
        scores = scores + " " +info[i] + "=" + data[i];
    }
    return scores;
  }


  // Populates the definition vcf lines ArrayList
  private void addDefinitionVcfLines() {
    for (String gene: m_definitionReader.getGenes()) {
      for (VariantLocus defPosition : m_definitionReader.getPositions(gene)) {
        // get each line of psuedo-vcf info from the definition file
        ArrayList<String> vcfLineDefinition = new ArrayList<>(Arrays.asList(
            ExtractPositions.getVcfLineFromDefinition(m_definitionReader, gene, defPosition, m_genomeBuild)));
        m_vcfDefinitionLine.put(m_definitionReader.getDefinitionFile(gene).getChromosome()  + "_" + defPosition.getVcfPosition(), vcfLineDefinition);
        }
      }
    }


  // Checks if position present in definition file
  private boolean checkPresent(String position) {
    boolean found = false;
    if (m_vcfDefinitionLine.containsKey(position)) {
      found = true;
    }
    return found;
  }


  // Lists vcf positions missing compared to the definition file
  private void checkMissingVcfLines() {
    System.out.println("\nThe following lines were missing from the vcf:");
    for (Map.Entry<String, ArrayList<String>> entry : m_vcfDefinitionLine.entrySet()) {
      if (!m_vcfMapFixed.containsKey( entry.getKey())) {
        System.out.print(String.join("\t", entry.getValue()));
        System.out.println("---->Position was missing from vcf");
        if (m_addMissing) {
          recordRepair(entry.getKey(), "Line was missing - added ref from def and presuming 0/0");
          entry.getValue().set(7, entry.getValue().get(7) + ";ADDED_MISSING;");
          m_vcfMapFixed.put(entry.getKey(), entry.getValue());
        }
      }
    }
  }

  // Writes vcf out
  private void writeVcf() {
    try (PrintWriter writer = new PrintWriter(String.valueOf(m_outputVcf), "UTF-8")) {
      for (String vcfHeader : m_vcfHeader) {
        writer.println(vcfHeader);
      }
      for (Map.Entry<String, ArrayList<String>> entry : m_vcfMapFixed.entrySet()) {
        String vcfLineOutput = String.join("\t", entry.getValue());
        writer.println(vcfLineOutput);
      }
      writer.flush();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  // Quick simple parsing of vcf:
  private void parseVcf()  {
    try {
      BufferedReader in = new BufferedReader(new FileReader(m_inputVcf.toString()));
      String line;
      while ((line = in.readLine()) != null) {
        if (!line.startsWith("#")) {
          //System.out.println(line);
          ArrayList<String> parts = new ArrayList<>(Arrays.asList(line.split("\t")));
          m_vcfMap.put(parts.get(0)+"_"+parts.get(1), parts);
        }
        else {
          m_vcfHeader.add(line);
        }
      }
      in.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  // HTML output table - cannot be too large
  private void vcfToHtml () {
    String resultsFile = m_inputVcf.toString()+".html";
    System.out.println("Printing html vcf for: " + resultsFile);
    StringBuilder builder = new StringBuilder();
    builder.append("<!DOCTYPE html>\n").append("<html>\n").append("<head>\n").append("<style>\n").append("table, th, td {\n")
        .append("    border: 1px solid black;\n").append("    border-collapse: collapse;\n").append("}\n").append("th, td {\n")
        .append("    padding: 5px;\n").append("    text-align: left;\n").append("}\n").append("table#t01 {\n").append("    width: 100%;    \n")
        .append("    background-color: #f1f1c1;\n").append("}\n").append("</style>\n")
        .append("<link href=\"https://netdna.bootstrapcdn.com/bootstrap/3.0.3/css/bootstrap.min.css\" rel=\"stylesheet\" /> ")
        .append("<link href=\"https://cdn.datatables.net/plug-ins/1.10.7/integration/bootstrap/3/dataTables.bootstrap.css\" rel=\"stylesheet\" /> ")
        .append("<link rel=\"stylesheet\" href=\"http://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.4.0/css/font-awesome.min.css\">")
        .append("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js\"></script> ")
        .append("<script src=\"https://cdn.datatables.net/1.10.7/js/jquery.dataTables.min.js\"></script> ")
        .append("<script src=\"https://cdn.datatables.net/plug-ins/1.10.7/integration/bootstrap/3/dataTables.bootstrap.js\"></script> ")
        .append("<script> $(document).ready(function() { ").append("$('#example').dataTable( { \"iDisplayLength\": 20 }); ").append("});")
        .append("</script>").append("<script>").append("$(document).ready( function() { ").append("$('table#example').show(); ")
        .append("$('i#spin').hide(); ").append("}); ").append("</script>").append("</head>\n").append("<body>").append("<div class=\"container\">");

    if (m_pharmcat) {
      builder.append(m_pharmCATResults);
    }

    builder.append("<H3>VCF Values (Please give about ten secs for long files to load <i id=\"spin\" class=\"fa fa-spinner fa-spin\" style=\"font-size:24px\"></i> ): </H3>");
    builder.append(" <table id=\"example\" class=\"table table-bordered\" cellspacing=\"0\" width=\"100%\" style=\"display:none\">\n" +
        "    <thead><tr> <th>chr</th> <th>Position</th> <th>rsid</th> <th>Ref</th> <th>Alt</th> <th>Genotype</th> <th>Star allele</th>" +
        "    <th>Info</th><th>Scores</th><th>Corrections</th> </tr>   </thead>\n<tbody>");

    for (Map.Entry<String, ArrayList<String>> entry : m_vcfMapFixed.entrySet()) {
      ArrayList<String> vcfLine = entry.getValue();

      String correction = "";
      if (m_vcfMapRepairs.containsKey(entry.getKey())) {
        correction = m_vcfMapRepairs.get(entry.getKey());
      }

      String positionString = "";
      if (m_vcfDefinitionLine.containsKey(entry.getKey())) {
        positionString = m_vcfDefinitionLine.get(entry.getKey()).get(7);
      }

      positionString = positionString.replace(",", ", ");
      String cellColor1= "";
      String cellColor2= "";
      String genotype = getGenotype(entry.getValue().get(8), entry.getValue().get(9));
      String info = entry.getValue().get(7).replace(";", "; ");
      String scores = getScores(entry.getValue().get(8), entry.getValue().get(9));

      if (genotype.contains("0/0") || genotype.contains("0|0")) {
        cellColor1 = "warning";
      } else if (genotype.contains("0/1") || genotype.contains("0|1")) {
        cellColor1= "success";
        cellColor2= "success";
      } else if (genotype.contains("1/1") || genotype.contains("1|1") ) {
        cellColor2= "success";
      }
      builder.append("\n<tr>" + "<td>").append(vcfLine.get(0)).append("</td> <td>").append(vcfLine.get(1)).append("</td> <td>")
          .append("<a href =\"https://www.pharmgkb.org/variant/").append(vcfLine.get(2)).append("\">").append(vcfLine.get(2))
          .append("</a>").append("</td> <td class=\"").append(cellColor1).append("\">").append(vcfLine.get(3).replace(",", ", "))
          .append("</td> <td class= \"").append(cellColor2).append("\">").append(vcfLine.get(4)).append("</td> <td>")
          .append(genotype).append("</td> <td>").append(positionString).append("</td> <td>").append(info).append("<td>")
          .append(scores).append("</td> <td>").append(correction).append("</td> </td> </tr>");
    }
    builder.append(" \n</tbody></table></div></body>");

    PrintWriter writer;
    try {
      writer = new PrintWriter(resultsFile, "UTF-8");
      writer.println(builder);
      writer.flush();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


}
