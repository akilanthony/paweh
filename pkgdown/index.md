<div class="paweh-home">

<section class="paweh-hero" aria-labelledby="paweh-title">
<div class="paweh-wordmark">paweh</div>
<h1 id="paweh-title">Power and Sample Size for<br class="paweh-hero-break"> Genetic Association Studies</h1>
<p class="paweh-lead">Plan genetic association studies before collecting data. Estimate power or the sample size you need under assumptions about genetic heterogeneity, genotype error, and phenotype misclassification.</p>
<div class="paweh-actions">
<a class="btn btn-primary" href="articles/paweh-01-getting-started.html">Get Started</a>
<a class="btn btn-outline-primary" href="#study-designs">Study Designs</a>
<a href="articles/paweh-08-methods-and-sources.html">Methods &amp; Sources</a>
<a href="reference/index.html">Function Reference</a>
</div>
</section>

<section class="paweh-section" aria-labelledby="planning-questions">
<h2 id="planning-questions">What can paweh calculate?</h2>
<p>Both calculations start with a proposed effect and a planned statistical test. In this <a href="articles/paweh-01-getting-started.html">worked example</a>, the two calculations use the same case-control assumptions and genotype chi-square test.</p>
<div class="paweh-example">
<div>
<h3>How much power does my study have?</h3>
<p class="paweh-result">98.52% power</p>
<p>With 500 cases and 500 controls, this is the chance the test detects the assumed association.</p>
</div>
<div>
<h3>How many participants do I need?</h3>
<p class="paweh-result">482 participants</p>
<p>For 80% power, recruit 241 cases and 241 controls. In a family study, plan for complete trios instead.</p>
</div>
</div>
</section>

<section class="paweh-section" aria-labelledby="study-designs">
<h2 id="study-designs">Choose your study design</h2>
<p>Start with the people you can recruit and the measurements you expect to collect.</p>
<ul class="paweh-designs">
<li>
<h3><a href="articles/paweh-02-case-control-study-design.html">Case-Control Studies</a></h3>
<p>Plan how many people with and without a condition to recruit for a comparison of genotypes.</p>
</li>
<li>
<h3><a href="articles/paweh-03-tdt-study-design.html">TDT / Family Studies</a></h3>
<p>Plan a study of affected children and both biological parents using informative allele transmissions.</p>
</li>
<li>
<h3><a href="articles/paweh-04-quantitative-trait-study-design.html">Quantitative Traits</a></h3>
<p>Plan a study of measurements such as height. Compare a full sample, selected extremes, or several traits.</p>
</li>
<li>
<h3><a href="articles/paweh-06-sequencing-study-design.html">Sequencing-Based Studies</a></h3>
<p>Compare how sequencing depth and read errors affect recruitment for case-control studies and trios.</p>
</li>
<li>
<h3><a href="articles/paweh-07-lrtae-study-design.html">Validation-Based / LRTae Studies</a></h3>
<p>Plan how many people to recruit when only some receive a reference diagnosis or a confirmatory genotype assay.</p>
</li>
<li>
<h3><a href="articles/paweh-09-hwe-diagnostics.html">Hardy–Weinberg Equilibrium</a></h3>
<p>Check genotype counts for departures from equilibrium and plan how many individuals you need to detect a specified departure.</p>
</li>
<li>
<h3><a href="articles/paweh-05-interactive-dashboard.html">Interactive Dashboard</a></h3>
<p>Explore case-control, trio, and quantitative-trait designs in a Shiny interface.</p>
</li>
</ul>
<p><a href="articles/index.html">Browse all articles</a></p>
</section>

<section class="paweh-section paweh-install" aria-labelledby="install">
<h2 id="install">Install and get started</h2>
<p>paweh is currently available from GitHub. The package has yet to be released on CRAN. Install the stable <code>main</code> branch with the command below.</p>
<pre><code>remotes::install_github("akilanthony/paweh", ref = "main")</code></pre>
<p>Install <code>remotes</code> first if needed. For the current development package and all nine articles, see the <code>develop</code> installation in <a href="articles/paweh-01-getting-started.html">Getting Started</a>.</p>
</section>

<section class="paweh-section" aria-labelledby="methodological-foundation">
<h2 id="methodological-foundation">Methodological foundation</h2>
<div class="paweh-book">
<figure>
<a href="https://doi.org/10.1007/978-3-030-61121-7"><img src="reference/figures/textbook-cover.jpg" width="620" height="934" alt="Springer cover of Heterogeneity in Statistical Genetics by Derek Gordon, Stephen J. Finch, and Wonkuk Kim" loading="lazy"></a>
<figcaption>Cover artwork by Springer. Used with permission.</figcaption>
</figure>
<div>
<p class="paweh-book-title"><cite>Heterogeneity in Statistical Genetics: How to Assess, Address, and Account for Mixtures in Association Studies</cite></p>
<p>Derek Gordon, Stephen J. Finch, and Wonkuk Kim. Springer, 2020.</p>
<p>The book explains how mixtures and measurement errors can change a genetic study’s power. Derek Gordon is a coauthor of both the book and paweh. The package also draws on other publications for individual methods.</p>
<p><a class="btn btn-outline-primary" href="articles/paweh-08-methods-and-sources.html">Read Methods &amp; Sources</a></p>
</div>
</div>
</section>

<section class="paweh-section paweh-scope" aria-labelledby="scope">
<h2 id="scope">Capabilities and limitations</h2>
<p>paweh supports prospective study design using specified effect sizes, sampling plans, and measurement assumptions. Its sequencing calculations use planned read depths and error rates. Use separate analysis tools to process raw sequencing files or run downstream GWAS association tests.</p>
<p><a href="authors.html#citation">Cite paweh</a></p>
</section>

</div>
