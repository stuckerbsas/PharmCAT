package org.pharmgkb.pharmcat.definition;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.junit.Test;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;

import static org.junit.Assert.assertEquals;


/**
 *  JUnit test for {@link ExtractIntervals}.
 * Created by lester on 10/14/16.
 */
public class ExtractIntervalsTest {

  // Quick test of extracting the intervals on a single definition file, making sure the count matches
  @Test
  public void testExtractIntervals() throws IOException {
    File file = new File(DefinitionReader.class.getResource("VKORC1_translation.json").getFile());
    DefinitionReader definitionReader = new DefinitionReader();
    Path path = Paths.get(file.getAbsolutePath());
    definitionReader.read(path);
    Path outputPath = Paths.get("org/pharmgkb/pharmcat/definition/intervals.list");
    ExtractIntervals extractIntervals = new ExtractIntervals(path,outputPath);
    StringBuilder extractIntervalsString = extractIntervals.getPositions(definitionReader);
    int definitionPositions = definitionReader.getPositions("VKORC1").length;
    int extractIntervalPositions = extractIntervalsString.toString().split("\n").length;
    assertEquals(definitionPositions, extractIntervalPositions);
  }


}
