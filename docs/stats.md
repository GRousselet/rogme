Statistical tests and measures of effect sizes in `rogme`
================
Guillaume A. Rousselet
2018-07-07

<table style="width:94%;">
<colgroup>
<col width="47%" />
<col width="47%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Function</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>hd(x,q=.5)</code></td>
<td align="left">Harrell-Davis quantile estimator</td>
</tr>
<tr class="even">
<td align="left"><code>hdseq(x, qseq=seq(0.1,0.9,0.1))</code></td>
<td align="left">Compute a sequence of quantiles</td>
</tr>
<tr class="odd">
<td align="left"><code>allpdiff(x,y)</code></td>
<td align="left">Calculate all pairwise differences between 2 vectors</td>
</tr>
<tr class="even">
<td align="left"><code>allpdiff_hdpbci(x,y,alpha=.05,q=.5,nboot=600)</code></td>
<td align="left">Compute percentile bootstrap confidence interval of the median of all pairwise differences</td>
</tr>
<tr class="odd">
<td align="left"><code>pb2gen(x,y,alpha=.05, nboot=2000, est=hd)</code></td>
<td align="left">Compare two independent groups using the percentile bootstrap - default estimator = <code>hd</code></td>
</tr>
<tr class="even">
<td align="left"><code>hdpbci(x,q=.5)</code></td>
<td align="left">Compute percentile bootstrap confidence interval of the qth quantile</td>
</tr>
<tr class="odd">
<td align="left"><code>quantiles_pbci(x,q=seq(.1,.9,.1),nboot=2000,alpha=0.05)</code></td>
<td align="left">Compute percentile bootstrap confidence intervals of a sequence of quantiles</td>
</tr>
<tr class="even">
<td align="left"><code>ks(x,y)</code></td>
<td align="left">Kolmogorov-Smirnov test statistic</td>
</tr>
<tr class="odd">
<td align="left"><code>cid(x,y)</code></td>
<td align="left">Cliff's delta test</td>
</tr>
<tr class="even">
<td align="left"><code>pxgta(x,a=0)</code></td>
<td align="left">Proportion of observations greater than a specified value</td>
</tr>
<tr class="odd">
<td align="left"><code>pxgty(x,y)</code></td>
<td align="left">Proportion of observations in x greater than observations in y</td>
</tr>
</tbody>
</table>
