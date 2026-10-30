// Builds the main-text and supplementary Word documents for the robustness analysis.
// usage: node make_docx.js   ->  manuscript/robust_analysis_main_text.docx, manuscript/robust_analysis_supplementary.docx
const fs = require("fs"), path = require("path");
const { Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType, LineRuleType } = require("docx");
const HERE = path.join(__dirname, ".."), FIG = path.join(HERE, "figures"), OUT = path.join(HERE, "manuscript");
const FONT = "Arial", SIZE = 22;                       // 11 pt
const TEXT_W_IN = 6.5;                                 // US Letter, 1-inch margins

function pngSize(file) { const b = fs.readFileSync(file); return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) }; }
// Text with inline markup: **bold**, *italic*, ^sup^, ~sub~
function runs(text) {
  const out = []; const re = /(\*\*[^*]+\*\*|\*[^*]+\*|\^[^^]+\^|~[^~]+~)/g; let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new TextRun({ text: text.slice(last, m.index) }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new TextRun({ text: t.slice(2, -2), bold: true }));
    else if (t.startsWith("*")) out.push(new TextRun({ text: t.slice(1, -1), italics: true }));
    else if (t.startsWith("^")) out.push(new TextRun({ text: t.slice(1, -1), superScript: true }));
    else out.push(new TextRun({ text: t.slice(1, -1), subScript: true }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new TextRun({ text: text.slice(last) }));
  return out;
}
const para = (text, opts = {}) => new Paragraph({ children: runs(text), spacing: { after: 160 }, alignment: AlignmentType.JUSTIFIED, ...opts });
const heading = (text) => new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text })], spacing: { before: 240, after: 160 } });
function figure(file) {
  const { w, h } = pngSize(path.join(FIG, file)); const W = TEXT_W_IN * 96;      // docx-js sizes images in pixels at 96 dpi
  return new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 120 }, keepNext: true,
    children: [new ImageRun({ type: "png", data: fs.readFileSync(path.join(FIG, file)), transformation: { width: W, height: Math.round(W * h / w) } })] });
}
function doc(children) {
  return new Document({
    styles: {
      default: { document: { run: { font: FONT, size: SIZE }, paragraph: { spacing: { line: 240, lineRule: LineRuleType.AUTO } } } },
      paragraphStyles: [{ id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { font: FONT, size: 26, bold: true }, paragraph: { spacing: { before: 240, after: 160 }, outlineLevel: 0 } }],
    },
    sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } }, children }],
  });
}

// ---------- shared caption pieces (style of the final supplementary, 2026-10-30) ----------
const MODEL = { standard: "phage–bacteria (Lotka–Volterra) model", cc: "phage–bacteria model with bacterial carrying capacity" };
const NPROC = { standard: 4, cc: 5 };
const SETDESC = {
  completed: "only simulations that remained within the physical bounds for the full observational window are shown",
  stopped: "only simulations that were stopped because a population left the physical bounds before the end of the observational window are shown, and both models are compared up to the stop time",
};
function robustCaption(label, title, model, set, n) {
  const P = NPROC[model], span = set === "stopped" ? "simulated" : "observational";
  return `**${label}. ${title}** Latin hypercube sampling (log-uniform) of the parameters, initial concentrations, and observational time of the ${MODEL[model]}; ${SETDESC[set]} (n = ${n} for α = 0.1). ` +
    `**a)** Histogram of the mean combined relative error (%) between the transient (Boolean) and full models at α = 0.1${set === "stopped" ? ", calculated over the simulated portion of each trajectory" : ""}; the bar labeled 0 counts runs with zero error (all processes relevant throughout), the dashed line marks the median, and the solid black line marks the observer tolerance (10%, corresponding to α = 0.1). ` +
    `**b)** Fraction of the ${span} time that each run spent with 0 to ${P} relevant processes (weight w ≥ 1) at α = 0.1; violins show the distribution across runs, gray points are a random sample of 200 runs, and the numbers on top (n) are the runs that spent any time with that number of relevant processes. ` +
    `**c)** Top, box plots of the mean combined relative error for α = 10^−5^ to 1 (box, interquartile range, IQR; line, median; whiskers, 1.5 × IQR; points, values beyond the whiskers), with the mean (◆) and the observer tolerance, 100·α % (dashed line); zero errors are drawn at the bottom of the axis, labeled 0. Bottom, mean ± standard deviation of the fraction of ${span} time with all ${P} processes relevant. ` +
    `Extreme outliers (Tukey's far-out test, above Q3 + 3 × IQR of log~10~ error) are excluded from the statistics and, when present, shown in ochre.`;
}

// ---------- main text ----------
const main = doc([
  heading("Results"),
  para(
    "**Robust analysis of the Boolean model.** To test whether simulating only the relevant processes reproduces the full dynamics across the parameter space, we sampled the model parameters, initial concentrations and observational times with Latin hypercube sampling and compared the Boolean (transient) model with the full model (Figure 6). " +
    "For a 10% observer tolerance (α = 0.1), the mean combined relative error was below the tolerance in 997 of 1,000 simulations, with a median of 0.68% and a 90th percentile of 2.2% (Figure 6a). " +
    "The number of relevant processes varied widely within and between trajectories (Figure 6b): on average, all four processes were relevant during 19% of the observational window, and only 9.5% of the trajectories (95% CI 7.8–11.5%) kept all processes relevant for the whole window, which represents 0.96% (0.82–1.13%) of the sampled space. " +
    "Across tolerances spanning five orders of magnitude (α = 10^−5^ to 1), the mean error remained below the observer tolerance in at least 99.7% of the simulations, and the error decreased as the tolerance decreased (Figure 6c). " +
    "As the tolerance increased, the fraction of time with all processes relevant decreased from 0.95 (α = 10^−5^) to 0.09 (α = 1), and the share of the sampled space requiring all processes for the whole window decreased from 9.0% to 0.4%. " +
    "The model with carrying capacity and the simulations that had to be stopped because a population left the physical bounds led to similar or better results (carrying capacity: median error 0.77%, 988 of 999 simulations within 10%; stopped simulations: median error 0.0035% without and 0.0065% with carrying capacity, more than 99.9% of simulations within 10%; Supplementary Figures S3, S4 and S5). " +
    "For α = 0.1, the carrying-capacity case showed 12 of 1,000 simulations above the tolerance (one of them an extreme outlier). These corresponded to long oscillatory trajectories in which the Boolean model accumulated a persistent delay with respect to the full model; the Boolean trajectories, relevant processes and process weights remained qualitatively in good agreement with the full model (Supplementary Figure S6)."),
  heading("Methods"),
  para(
    "**Robust analysis.** The full and Boolean phage–bacteria models, with and without carrying capacity, were taken unchanged from the project notebook and integrated with *solve_ivp* (LSODA; relative and absolute tolerances 10^−13^ and 10^−16^) on 10,000 evenly spaced time points over the observational time, τ. In the Boolean model, a process is active when its weight is w ≥ 1. " +
    "Model parameters (r, a, c, m) and τ were sampled within the ranges in Table 1, initial concentrations within 10^3^–10^12^ bacteria ml^−1^ and 10^3^–10^14^ phages ml^−1^, and the carrying capacity within K = 10^4^–10^13^ cells ml^−1^ with B~0~ ≤ K, using log-uniform Latin hypercube sampling (scipy.stats.qmc; seed 2026). " +
    "Each sample was simulated with the same parameters at α = 10^−5^, 10^−4^, …, 1. A simulation was stopped when either population, in either model, fell below 10^−3^ ml^−1^ or exceeded the physical limits of 10^15^ bacteria ml^−1^ or 10^18^ phages ml^−1^ (the volumes of the smallest bacterium, about 10^−3^ µm^3^, and phage, about 10^−6^ µm^3^, packed into 1 ml); both models were then compared up to the earliest stop time. " +
    "The sample size was doubled (2,000, 4,000, …) until every α had at least 1,000 simulations covering the full observational window (16,000 samples without and 8,000 with carrying capacity), and 1,000 of these were selected at random for each α; stopped simulations with at least two compared time points were analyzed separately. " +
    "The accuracy of the Boolean model was measured as the mean combined relative error, the time average of ½(|B~full~ − B~Boolean~|/B~full~ + |P~full~ − P~Boolean~|/P~full~) × 100. For each simulation, we also recorded the fraction of time spent with 0 to P relevant processes (P = 4 without and 5 with carrying capacity) and whether all processes were relevant throughout. " +
    "Extreme outliers were identified with Tukey's far-out test on log~10~ error (values above Q3 + 3 × IQR) and excluded from the summary statistics, and 95% confidence intervals for proportions used the Wilson score interval. " +
    "The code and data are available in the project repository (results/alpha-robustness-v2)."),
  heading("Figure"),
  figure("standard_completed_figure.png"),
  para(robustCaption("Figure 6", "Robust analysis of the Boolean model.", "standard", "completed", "1,000")),
]);

// ---------- supplementary ----------
const supp = doc([
  heading("Figures"),
  figure("cc_completed_figure.png"),
  para(robustCaption("Figure S3", "Robust analysis of the Boolean model with carrying capacity.", "cc", "completed", "1,000")),
  figure("standard_stopped_figure.png"),
  para(robustCaption("Figure S4", "Robust analysis of the Boolean model for stopped simulations.", "standard", "stopped", "6,290")),
  figure("cc_stopped_figure.png"),
  para(robustCaption("Figure S5", "Robust analysis of the Boolean model with carrying capacity for stopped simulations.", "cc", "stopped", "3,513")),
  figure("figure_S6_cc_outlier_run178.png"),
  para("**Figure S6. Example of the extreme outlier with carrying capacity: persistent delay.** Simulation 178 of the robust analysis of the model with carrying capacity at α = 0.1 (mean combined relative error 8.6 × 10^5^%; extreme outlier). Blue, bacteria and bacterial processes; red, phage and phage processes; dashed black, transient (Boolean) model. " +
    "**a)** Bacteria and phage concentrations over the full observational window for the full (solid) and transient (dashed) models; the shaded interval marks the last cycles shown in b–d. " +
    "**b)** Concentrations over the last cycles, showing the delay accumulated by the transient model. " +
    "**c)** Number of relevant processes (w ≥ 1) in the full (solid gray) and transient (dashed) models. " +
    "**d)** Weight of each process in the full (colored) and transient (dashed) models; the gray band marks irrelevant weights (w < 1), and the solid line marks the threshold w = 1."),
]);

fs.mkdirSync(OUT, { recursive: true });
Packer.toBuffer(main).then(b => fs.writeFileSync(path.join(OUT, "robust_analysis_main_text.docx"), b));
Packer.toBuffer(supp).then(b => fs.writeFileSync(path.join(OUT, "robust_analysis_supplementary.docx"), b));
