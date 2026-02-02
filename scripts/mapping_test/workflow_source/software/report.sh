#!/bin/bash

sampleName="${1}"
preFlagstatCollapse="${2}"
preStatsCollapse="${3}"
preDepthDistributionCollapse="${4}"
preDepthAlongReferenceCollapse="${5}"
postFlagstatCollapse="${6}"
postStatsCollapse="${7}"
postDepthDistributionCollapse="${8}"
postDepthAlongReferenceCollapse="${9}"
preFlagstatDiscard="${10}"
preStatsDiscard="${11}"
preDepthDistributionDiscard="${12}"
preDepthAlongReferenceDiscard="${13}"
postFlagstatDiscard="${14}"
postStatsDiscard="${15}"
postDepthDistributionDiscard="${16}"
postDepthAlongReferenceDiscard="${17}"
preFlagstatIgnore="${18}"
preStatsIgnore="${19}"
preDepthDistributionIgnore="${20}"
preDepthAlongReferenceIgnore="${21}"
postFlagstatIgnore="${22}"
postStatsIgnore="${23}"
postDepthDistributionIgnore="${24}"
postDepthAlongReferenceIgnore="${25}"
outputDir="${26}"

cat << EOF > "$outputDir"/report.css
table, th, td {
	border: 1px solid black;
}
EOF

cat << EOF > "$outputDir"/"$sampleName".report.html
<!DOCTYPE html>
<html>
	<head>
		<link rel="stylesheet" href="report.css" type="css"/>
		<title>Data Report</title>
	</head>
	<body>
		<h1>Data Report</h1>
		<h2>$sampleName</h2>
		<h3>Summary Statistics</h3>
		<h4>Pre-filtering</h4>
		<table>
			<tr>
				<th class="column1">

				</th>
				<th class="column2">
					Collapse
				</th>
				<th class="column3">
					Discard
				</th>
				<th class="column4">
					Ignore
				</th>
			</tr>
			<tr>
				<td class="column1">
					Total number of reads
				</td>
				<td class="column2">
					$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Primary reads
				</td>
				<td class="column2">
					$(grep "primary$" "$preFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "primary$" "$preFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "primary$" "$preFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Secondary reads
				</td>
				<td class="column2">
					$(grep "secondary" "$preFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "secondary" "$preFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "secondary" "$preFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Supplementary reads
				</td>
				<td class="column2">
					$(grep "supplementary" "$preFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "supplementary" "$preFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "supplementary" "$preFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Duplicate reads
				</td>
				<td class="column2">
					$(grep "reads duplicated:" "$preStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "reads duplicated:" "$preStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "reads duplicated:" "$preStatsIgnore" | cut -f 3)
				</td>
			</tr>
			<tr>
				<td class="column1">
					Mapped reads
				</td>
				<td class="column2">
					$(grep "reads mapped:" "$preStatsCollapse" | cut -f 3) / $(awk -v mapped="$(grep "reads mapped:" "$preStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print mapped / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads mapped:" "$preStatsDiscard" | cut -f 3) / $(awk -v mapped="$(grep "reads mapped:" "$preStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print mapped / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads mapped:" "$preStatsIgnore" | cut -f 3) / $(awk -v mapped="$(grep "reads mapped:" "$preStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print mapped / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Primary mapped reads
				</td>
				<td class="column2">
					$(grep "primary mapped" "$preFlagstatCollapse" | sed 's/ \+.*//') / $(awk -v primaryMapped="$(grep "primary mapped" "$preFlagstatCollapse" | sed 's/ \+.*//')" -v total="$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print primaryMapped / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "primary mapped" "$preFlagstatDiscard" | sed 's/ \+.*//') / $(awk -v primaryMapped="$(grep "primary mapped" "$preFlagstatDiscard" | sed 's/ \+.*//')" -v total="$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print primaryMapped / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "primary mapped" "$preFlagstatIgnore" | sed 's/ \+.*//') / $(awk -v primaryMapped="$(grep "primary mapped" "$preFlagstatIgnore" | sed 's/ \+.*//')" -v total="$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print primaryMapped / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Unmapped reads
				</td>
				<td class="column2">
					$(grep "reads unmapped:" "$preStatsCollapse" | cut -f 3) / $(awk -v unmapped="$(grep "reads unmapped:" "$preStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print unmapped / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads unmapped:" "$preStatsDiscard" | cut -f 3) / $(awk -v unmapped="$(grep "reads unmapped:" "$preStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print unmapped / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads unmapped:" "$preStatsIgnore" | cut -f 3) / $(awk -v unmapped="$(grep "reads unmapped:" "$preStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print unmapped / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Paired reads
				</td>
				<td class="column2">
					$(grep "reads paired:" "$preStatsCollapse" | cut -f 3) / $(awk -v paired="$(grep "reads paired:" "$preStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print paired / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads paired:" "$preStatsDiscard" | cut -f 3) / $(awk -v paired="$(grep "reads paired:" "$preStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print paired / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads paired:" "$preStatsIgnore" | cut -f 3) / $(awk -v paired="$(grep "reads paired:" "$preStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print paired / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Properly paired reads
				</td>
				<td class="column2">
					$(grep "reads properly paired:" "$preStatsCollapse" | cut -f 3) / $(awk -v properlyPaired="$(grep "reads properly paired:" "$preStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print properlyPaired / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads properly paired:" "$preStatsDiscard" | cut -f 3) / $(awk -v properlyPaired="$(grep "reads properly paired:" "$preStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print properlyPaired / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads properly paired:" "$preStatsIgnore" | cut -f 3) / $(awk -v properlyPaired="$(grep "reads properly paired:" "$preStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print properlyPaired / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Singletons
				</td>
				<td class="column2">
					$(grep "singletons" "$preFlagstatCollapse" | sed 's/ \+.*//') / $(awk -v single="$(grep "singletons" "$preFlagstatCollapse" | sed 's/ \+.*//')" -v total="$(grep "total" "$preFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print single / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "singletons" "$preFlagstatDiscard" | sed 's/ \+.*//') / $(awk -v single="$(grep "singletons" "$preFlagstatDiscard" | sed 's/ \+.*//')" -v total="$(grep "total" "$preFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print single / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "singletons" "$preFlagstatIgnore" | sed 's/ \+.*//') / $(awk -v single="$(grep "singletons" "$preFlagstatIgnore" | sed 's/ \+.*//')" -v total="$(grep "total" "$preFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print single / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Average read length
				</td>
				<td class="column2">
					$(grep "average length:" "$preStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "average length:" "$preStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "average length:" "$preStatsIgnore" | cut -f 3)
				</td>
			</tr>
			<tr>
				<td class="column1">
					Average insert size
				</td>
				<td class="column2">
					$(grep "insert size average:" "$preStatsCollapse" | cut -f 3) +- $(grep "insert size standard deviation:" "$preStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "insert size average:" "$preStatsDiscard" | cut -f 3) +- $(grep "insert size standard deviation:" "$preStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "insert size average:" "$preStatsIgnore" | cut -f 3) +- $(grep "insert size standard deviation:" "$preStatsIgnore" | cut -f 3)
				</td>
			</tr>
			<tr>
				<td class="column1">
					Average mapping quality
				</td>
				<td class="column2">
					$(grep "average quality:" "$preStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "average quality:" "$preStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "average quality:" "$preStatsIgnore" | cut -f 3)
				</td>
			</tr>
		</table>
		<h4>Post-filtering</h4>
		<table>
			<tr>
				<th class="column1">

				</th>
				<th class="column2">
					Collapse
				</th>
				<th class="column3">
					Discard
				</th>
				<th class="column4">
					Ignore
				</th>
			</tr>
			<tr>
				<td class="column1">
					Total number of reads
				</td>
				<td class="column2">
					$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Primary reads
				</td>
				<td class="column2">
					$(grep "primary$" "$postFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "primary$" "$postFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "primary$" "$postFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Secondary reads
				</td>
				<td class="column2">
					$(grep "secondary" "$postFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "secondary" "$postFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "secondary" "$postFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Supplementary reads
				</td>
				<td class="column2">
					$(grep "supplementary" "$postFlagstatCollapse" | sed 's/ \+.*//')
				</td>
				<td class="column3">
					$(grep "supplementary" "$postFlagstatDiscard" | sed 's/ \+.*//')
				</td>
				<td class="column4">
					$(grep "supplementary" "$postFlagstatIgnore" | sed 's/ \+.*//')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Duplicate reads
				</td>
				<td class="column2">
					$(grep "reads duplicated:" "$postStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "reads duplicated:" "$postStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "reads duplicated:" "$postStatsIgnore" | cut -f 3)
				</td>
			</tr>
			<tr>
				<td class="column1">
					Mapped reads
				</td>
				<td class="column2">
					$(grep "reads mapped:" "$postStatsCollapse" | cut -f 3) / $(awk -v mapped="$(grep "reads mapped:" "$postStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print mapped / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads mapped:" "$postStatsDiscard" | cut -f 3) / $(awk -v mapped="$(grep "reads mapped:" "$postStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print mapped / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads mapped:" "$postStatsIgnore" | cut -f 3) / $(awk -v mapped="$(grep "reads mapped:" "$postStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print mapped / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Primary mapped reads
				</td>
				<td class="column2">
					$(grep "primary mapped" "$postFlagstatCollapse" | sed 's/ \+.*//') / $(awk -v primaryMapped="$(grep "primary mapped" "$postFlagstatCollapse" | sed 's/ \+.*//')" -v total="$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print primaryMapped / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "primary mapped" "$postFlagstatDiscard" | sed 's/ \+.*//') / $(awk -v primaryMapped="$(grep "primary mapped" "$postFlagstatDiscard" | sed 's/ \+.*//')" -v total="$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print primaryMapped / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "primary mapped" "$postFlagstatIgnore" | sed 's/ \+.*//') / $(awk -v primaryMapped="$(grep "primary mapped" "$postFlagstatIgnore" | sed 's/ \+.*//')" -v total="$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print primaryMapped / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Unmapped reads
				</td>
				<td class="column2">
					$(grep "reads unmapped:" "$postStatsCollapse" | cut -f 3) / $(awk -v unmapped="$(grep "reads unmapped:" "$postStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print unmapped / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads unmapped:" "$postStatsDiscard" | cut -f 3) / $(awk -v unmapped="$(grep "reads unmapped:" "$postStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print unmapped / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads unmapped:" "$postStatsIgnore" | cut -f 3) / $(awk -v unmapped="$(grep "reads unmapped:" "$postStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print unmapped / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Paired reads
				</td>
				<td class="column2">
					$(grep "reads paired:" "$postStatsCollapse" | cut -f 3) / $(awk -v paired="$(grep "reads paired:" "$postStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print paired / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads paired:" "$postStatsDiscard" | cut -f 3) / $(awk -v paired="$(grep "reads paired:" "$postStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print paired / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads paired:" "$postStatsIgnore" | cut -f 3) / $(awk -v paired="$(grep "reads paired:" "$postStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print paired / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Properly paired reads
				</td>
				<td class="column2">
					$(grep "reads properly paired:" "$postStatsCollapse" | cut -f 3) / $(awk -v properlyPaired="$(grep "reads properly paired:" "$postStatsCollapse" | cut -f 3)" -v total="$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print properlyPaired / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "reads properly paired:" "$postStatsDiscard" | cut -f 3) / $(awk -v properlyPaired="$(grep "reads properly paired:" "$postStatsDiscard" | cut -f 3)" -v total="$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print properlyPaired / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "reads properly paired:" "$postStatsIgnore" | cut -f 3) / $(awk -v properlyPaired="$(grep "reads properly paired:" "$postStatsIgnore" | cut -f 3)" -v total="$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print properlyPaired / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Singletons
				</td>
				<td class="column2">
					$(grep "singletons" "$postFlagstatCollapse" | sed 's/ \+.*//') / $(awk -v single="$(grep "singletons" "$postFlagstatCollapse" | sed 's/ \+.*//')" -v total="$(grep "total" "$postFlagstatCollapse" | sed 's/ \+.*//')" 'BEGIN{print single / total * 100 "%"}')
				</td>
				<td class="column3">
					$(grep "singletons" "$postFlagstatDiscard" | sed 's/ \+.*//') / $(awk -v single="$(grep "singletons" "$postFlagstatDiscard" | sed 's/ \+.*//')" -v total="$(grep "total" "$postFlagstatDiscard" | sed 's/ \+.*//')" 'BEGIN{print single / total * 100 "%"}')
				</td>
				<td class="column4">
					$(grep "singletons" "$postFlagstatIgnore" | sed 's/ \+.*//') / $(awk -v single="$(grep "singletons" "$postFlagstatIgnore" | sed 's/ \+.*//')" -v total="$(grep "total" "$postFlagstatIgnore" | sed 's/ \+.*//')" 'BEGIN{print single / total * 100 "%"}')
				</td>
			</tr>
			<tr>
				<td class="column1">
					Average read length
				</td>
				<td class="column2">
					$(grep "average length:" "$postStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "average length:" "$postStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "average length:" "$postStatsIgnore" | cut -f 3)
				</td>
			</tr>
			<tr>
				<td class="column1">
					Average insert size
				</td>
				<td class="column2">
					$(grep "insert size average:" "$postStatsCollapse" | cut -f 3) +- $(grep "insert size standard deviation:" "$postStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "insert size average:" "$postStatsDiscard" | cut -f 3) +- $(grep "insert size standard deviation:" "$postStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "insert size average:" "$postStatsIgnore" | cut -f 3) +- $(grep "insert size standard deviation:" "$postStatsIgnore" | cut -f 3)
				</td>
			</tr>
			<tr>
				<td class="column1">
					Average mapping quality
				</td>
				<td class="column2">
					$(grep "average quality:" "$postStatsCollapse" | cut -f 3)
				</td>
				<td class="column3">
					$(grep "average quality:" "$postStatsDiscard" | cut -f 3)
				</td>
				<td class="column4">
					$(grep "average quality:" "$postStatsIgnore" | cut -f 3)
				</td>
			</tr>
		</table>
		<h3>Depth Distribution</h3>
		<h4>Pre-filtering</4>
		<h5>Collapse</5>
		<img class="graph" src="$preDepthDistributionCollapse" width="600"/>
		<h5>Discard</5>
		<img class="graph" src="$preDepthDistributionDiscard" width="600"/>
		<h5>Ignore</5>
		<img class="graph" src="$preDepthDistributionIgnore" width="600"/>
		<h4>Post-filtering</4>
		<h5>Collapse</h5>
		<img class="graph" src="$postDepthDistributionCollapse" width="600"/>
		<h5>Discard</h5>
		<img class="graph" src="$postDepthDistributionDiscard" width="600"/>
		<h5>Ignore</h5>
		<img class="graph" src="$postDepthDistributionIgnore" width="600"/>
		<h3>Depth Along Reference</h3>
		<h4>Pre-filtering</4>
		<h5>Collapse</h5>
		<img class="graph" src="$preDepthAlongReferenceCollapse" width="600"/>
		<h5>Discard</h5>
		<img class="graph" src="$preDepthAlongReferenceDiscard" width="600"/>
		<h5>Ignore</h5>
		<img class="graph" src="$preDepthAlongReferenceIgnore" width="600"/>
		<h4>Post-filtering</4>
		<h5>Collapse</h5>
		<img class="graph" src="$postDepthAlongReferenceCollapse" width="600"/>
		<h5>Discard</h5>
		<img class="graph" src="$postDepthAlongReferenceDiscard" width="600"/>
		<h5>Ignore</h5>
		<img class="graph" src="$postDepthAlongReferenceIgnore" width="600"/>
	</body>
</html>
EOF