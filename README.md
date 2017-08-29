## Regional variation in DNA gain and loss

__A collection of custom scripts designed to calculate regional levels of DNA gain and loss between two species.__

##### getGapsNet
Takes UCSC net file as input and extracts gaps.
Returns reference and query start and end coordinates for each gap as output.
	
##### getFillsNet
Takes UCSC net file as input and extracts fills.
Returns reference and query start and end coordinates for each fill as output.

getFillsNet is also used to extract fills between the reference genome and each outgroup species.
Outgroup species fills are further processed downstream  using bedtools to create our 'ancestral elements'.

##### formatData.R
Takes several inputs for both a reference and query genome:

+ Gaps
+ Fills
+ Ancestral elements

After specifying the reference (ref) and query (que) genome, relevant genome annotations are downloaded directly from the UCSC database.
These include `chromInfo.txt` which contains chromosome sizes for a particular assembly and `gaps.txt` which contains coordinates for assembly gaps.
Input data is reformatted into GRanges objects where each ref interval for gaps and fills is paired with que information and vice versa.
Non-reciprocal best hit (RBH) for gaps (alignment gaps, not assembly gaps) and fills are removed.
Resulting output consists of several R objects compressed into a single file named `ref.que.netData.RData`.
The compressed R object output includes:

+ fills (RBH only)
+ gaps (RBH only)
+ ancestral elements
+ assembly gaps
+ chromosome lengths
+ non-RBH gaps (removed from 'gaps' data)
+ non-RBH gaps (removed from 'fills' data)
+ fill gaps (all gaps between fills, includes those outside of nets)


##### sortRepeatFamilies.R

Repeat Masker annotations for ref and que along with output from formatData.R are received as input.
Ancestral elements and fills are used to identify repeat families with lineage specific activity.
Recent transposon families to be used to identify DNA gain is returned as output.


##### syntheticGenomeBuilder.R
Output from formatData.R is used as input.
There is also an option to include recent transposon output from sortRepeatFamilies.R too.
Gaps in ref and que are annotated as gain and loss using one of two methods depending on whether sortRepeatFamilies.R is provided.

The first method for gap annotation uses ancestral elements.
For example, que gaps overlapping an ancestral element are identified as que loss.
Meanwhile que gaps not overlapping an ancestral element are identified as ref gain.
The second method for gap annotation uses recent transposons.
In this case, que gaps overlapping a recent transposon are identified as ref gain.
Meanwhile que gaps not overlapping a recent transposon are identified as que loss.
Each method is complementary and relies on the correct identification of ancestral elements and recent transposons.

Next, All annotated gaps are placed within the ref genomic background.
Positions in ref downstream of recently placed que gap annotations are shifted downstream relative to the size of the placed gap.
Finally, A stretched ref genome containing gap annotations from both ref and que is returned as output.

##### synthChrom.R
Stretched genomes produced as output from syntheticGenomeBuilder.R are received as input.
Stretched genomes are broken up into bins of 200kb and the ref gain, ref loss, que gain and que loss is tallied for each bin.
A binned genome is produced as output.

##### methodComparisonHotspotID.R
Recent transposon based stretched genome and ancestral element based stretched genome are received as input.
Ref gain, ref loss, que gain and que loss hotspots are identified for both stretched genomes.
Hotspot bins are identified using the G*i statistics that takes neighbouring bins into account.
Genomic coordinates for hotspots bins are returned as output.

##### hotspotFeatureEnrichment.R
Hotspot bins and binned distribution of various genomic features are received as input.
10000 permutations are used to calculate a background distribution for each feature.
This is used to calculate level of enrichment for each feature in hotspots. 
Enrichment levels are returned as output in the form of Z scores.

##### GOtermAnalysis.R
Hotspots were received as input.
Genes were extracted and topGO was used to identify enriched GO terms.
Enrichment P-values were returned as output.



##### netDataFunctions.R
Functions frequently used by each script.