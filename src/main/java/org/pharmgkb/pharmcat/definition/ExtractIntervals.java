package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import javax.annotation.Nonnull;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;


/**
 * This extracts the positions from the definition files and reformats them as a GATK interval file, with a five nucleotide
 * padding either side.
 * As a command line it takes three arguments:
 * -d  CPIC definitions directory
 * -o output
 * -g genome build to use for DAS
 *
 * Created by lester on 08/02/16.
 */

public class ExtractIntervals {
  private static Path ps_definitionDir;
  private static Path ps_outputVcf;


  // Default constructor
  public ExtractIntervals(Path definitionDir, Path outputVcf) {
    ps_definitionDir = definitionDir;
    ps_outputVcf = outputVcf;
  }


  // Main method for command line use
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("d", "definition-dir", "directory of allele definition files", true, "d")
          .addOption("o", "output-file", "output vcf file", true, "o");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path definitionDir = cliHelper.getValidDirectory("d", false);
      Path outputVcf= cliHelper.getValidFile("o", false);


      ExtractIntervals extractPositions = new ExtractIntervals(definitionDir, outputVcf);
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
      }
      StringBuilder intervalString = getPositions(definitionReader);
      try (PrintWriter writer = new PrintWriter(String.valueOf(ps_outputVcf), "UTF-8")) {  // close PrintWriter with try
        writer.print(intervalString);
        writer.flush();
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }

  }


  // Build up string
  public  StringBuilder getPositions(@Nonnull DefinitionReader definitionReader) throws IOException {
    StringBuilder builder = new StringBuilder();
    for (String gene : definitionReader.getGenes()) {  // for each definition file
      for (VariantLocus variantLocus :definitionReader.getPositions(gene)) {
          int start = variantLocus.getVcfPosition()-5;
          int end = variantLocus.getVcfPosition()+5;
          String interval = definitionReader.getDefinitionFile(gene).getChromosome() +
              ":" +
              start + "-" + end;
          System.out.println("Gene:" + gene + " - " + interval);
          builder.append(interval).append("\n");
      }
    }
    return  builder;
  }
}

