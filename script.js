// ---- shared helpers ----
const MW_PER_BP = 650; // g/mol per bp (average, dsDNA)

// ng = pmol * bp * 650 / 1000
function ngFromPmol(pmol, bp) {
  return (pmol * bp * MW_PER_BP) / 1000;
}

// pmol = ng * 1000 / (bp * 650)
function pmolFromNg(ng, bp) {
  return (ng * 1000) / (bp * MW_PER_BP);
}

function fmt(n, digits = 2) {
  if (!isFinite(n)) return "—";
  return n.toLocaleString(undefined, { maximumFractionDigits: digits });
}

// convert a length input + unit ("bp"/"kb") to bp
function bpFromLength(value, unit) {
  const n = parseFloat(value);
  if (!n || n <= 0) return NaN;
  return unit === "kb" ? n * 1000 : n;
}

// ---- master mix solver ----
// fragments: [{ label, targetNg, stockConc }]
// Picks the largest pipetting-volume floor (between fallbackMin and idealMin) that still
// lets every fragment fit inside dnaVol. Fragments whose raw (undiluted) volume already
// clears the floor are pipetted directly; fragments below the floor are diluted up to it.
// If even fallbackMin doesn't fit, the DNA/water (and MM) volume is scaled up instead.
function solveMasterMix(fragments, dnaVol, idealMin, fallbackMin) {
  const info = fragments.map((f) => ({ ...f, rawVol: f.targetNg / f.stockConc }));
  const totalAt = (m) => info.reduce((sum, f) => sum + Math.max(f.rawVol, m), 0);

  let floor, effectiveDnaVol, status;
  const idealTotal = totalAt(idealMin);

  if (idealTotal <= dnaVol) {
    floor = idealMin;
    effectiveDnaVol = dnaVol;
    status = "ideal";
  } else {
    const fallbackTotal = totalAt(fallbackMin);
    if (fallbackTotal <= dnaVol) {
      let lo = fallbackMin, hi = idealMin;
      for (let i = 0; i < 40; i++) {
        const mid = (lo + hi) / 2;
        if (totalAt(mid) <= dnaVol) lo = mid; else hi = mid;
      }
      floor = lo;
      effectiveDnaVol = dnaVol;
      status = "reduced";
    } else {
      floor = fallbackMin;
      effectiveDnaVol = fallbackTotal;
      status = "scaled";
    }
  }

  const rows = info.map((f) => {
    const vol = Math.max(f.rawVol, floor);
    const diluted = f.rawVol < floor - 1e-9;
    const dilutedConc = diluted ? f.targetNg / vol : f.stockConc;
    const factor = diluted ? f.stockConc / dilutedConc : 1;
    return { ...f, vol, diluted, dilutedConc, factor };
  });

  const fragTotal = rows.reduce((sum, r) => sum + r.vol, 0);
  const water = Math.max(effectiveDnaVol - fragTotal, 0);
  const mmVol = effectiveDnaVol;

  return { rows, water, mmVol, dnaVol: effectiveDnaVol, totalReaction: effectiveDnaVol + mmVol, status, floor };
}

function masterMixStatusText(result, idealMin, fallbackMin) {
  if (result.status === "ideal") {
    return `All fragments pipetted at ≥ ${fmt(idealMin)} µL — no dilution adjustments needed.`;
  }
  if (result.status === "reduced") {
    return `Some fragment volumes reduced toward ${fmt(fallbackMin)} µL to fit the ${fmt(result.dnaVol)} µL DNA/water budget.`;
  }
  return `DNA/water volume increased to ${fmt(result.dnaVol)} µL (2× HiFi MM matched) to fit all fragments at the ${fmt(fallbackMin)} µL minimum pipetting volume.`;
}

function renderMmVolCell(cell, row) {
  if (!row.diluted) {
    cell.innerHTML = `${fmt(row.vol)} µL`;
    return;
  }
  cell.innerHTML =
    `${fmt(row.vol)} µL<span class="mm-note">${fmt(row.dilutedConc)} ng/µL (dilute 1:${fmt(row.factor)})</span>`;
}

// ---- top-level tool nav (placeholder for future calculators) ----
const printToolName = document.getElementById("print-tool-name");

document.querySelectorAll(".tool-nav-btn[data-tool]").forEach((btn) => {
  btn.addEventListener("click", () => {
    document.querySelectorAll(".tool-nav-btn[data-tool]").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".tool-panel").forEach((p) => p.classList.remove("active"));
    btn.classList.add("active");
    const panel = document.getElementById(`tool-${btn.dataset.tool}`);
    panel.classList.add("active");
    const printBtn = panel.querySelector(".print-btn");
    if (printBtn) printToolName.textContent = printBtn.dataset.printTitle;
  });
});

// ---- print / export ----
const printDate = document.getElementById("print-date");
function stampPrintDate() {
  printDate.textContent = `Printed ${new Date().toLocaleString(undefined, {
    dateStyle: "medium",
    timeStyle: "short",
  })}`;
}

document.querySelectorAll(".print-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    stampPrintDate();
    window.print();
  });
});

// also stamp the date if the user prints via the browser's own menu/Cmd+P
window.addEventListener("beforeprint", stampPrintDate);

// ---- mode switch (2-fragment vs multi-fragment) ----
// scoped to #hifi-mode-switch specifically — other "I know X" toggles also use .mode-btn
// for shared styling, and a plain document-wide `.mode-btn` selector here would collide
// with them (and with each other) since they don't share this switch's data-mode/mode-panel wiring.
document.getElementById("hifi-mode-switch").querySelectorAll(".mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    document.getElementById("hifi-mode-switch").querySelectorAll(".mode-btn").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".mode-panel").forEach((p) => p.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById(`mode-${btn.dataset.mode}`).classList.add("active");
  });
});

// ================= 2-FRAGMENT MODE =================
const twoVectorLen = document.getElementById("two-vector-len");
const twoVectorUnit = document.getElementById("two-vector-unit");
const twoInsertLen = document.getElementById("two-insert-len");
const twoInsertUnit = document.getElementById("two-insert-unit");
const twoRatio = document.getElementById("two-ratio");
const twoKnownNg = document.getElementById("two-known-ng");
const twoKnownTitle = document.getElementById("two-known-title");
const twoKnownLabel = document.getElementById("two-known-label");
const twoResultLabel = document.getElementById("two-result-label");
const twoResultNg = document.getElementById("two-result-ng");
const twoResultDetail = document.getElementById("two-result-detail");
const twoFormula = document.getElementById("two-formula");
const twoKnowSwitch = document.getElementById("two-know-switch");

let twoKnow = "vector"; // "vector" = have vector amount, solve for insert. "insert" = reverse.

twoKnowSwitch.querySelectorAll(".mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    twoKnowSwitch.querySelectorAll(".mode-btn").forEach((b) => b.classList.remove("active"));
    btn.classList.add("active");
    twoKnow = btn.dataset.know;

    if (twoKnow === "vector") {
      twoKnownTitle.textContent = "Vector amount";
      twoKnownLabel.firstChild.textContent = "Amount you have (ng)";
      twoResultLabel.textContent = "Insert needed";
      twoFormula.textContent = "insert (ng) = ratio × vector (ng) × (insert bp / vector bp)";
    } else {
      twoKnownTitle.textContent = "Insert amount";
      twoKnownLabel.firstChild.textContent = "Amount you have (ng)";
      twoResultLabel.textContent = "Vector needed";
      twoFormula.textContent = "vector (ng) = (insert (ng) / ratio) × (vector bp / insert bp)";
    }
    calcTwoFragment();
  });
});

// returns { vectorBp, insertBp, vectorNg, insertNg, ratio } or null if inputs are incomplete
function getTwoFragmentTargets() {
  const vectorBp = bpFromLength(twoVectorLen.value, twoVectorUnit.value);
  const insertBp = bpFromLength(twoInsertLen.value, twoInsertUnit.value);
  const ratio = parseFloat(twoRatio.value);
  const knownNg = parseFloat(twoKnownNg.value);

  if (!vectorBp || !insertBp || !knownNg || isNaN(ratio) || ratio <= 0 || knownNg <= 0) {
    return null;
  }

  if (twoKnow === "vector") {
    const vectorPmol = pmolFromNg(knownNg, vectorBp);
    const insertPmol = ratio * vectorPmol;
    const insertNg = ngFromPmol(insertPmol, insertBp);
    return { vectorBp, insertBp, vectorNg: knownNg, insertNg, ratio, vectorPmol, insertPmol };
  }

  const insertPmol = pmolFromNg(knownNg, insertBp);
  const vectorPmol = insertPmol / ratio;
  const vectorNg = ngFromPmol(vectorPmol, vectorBp);
  return { vectorBp, insertBp, vectorNg, insertNg: knownNg, ratio, vectorPmol, insertPmol };
}

function calcTwoFragment() {
  const t = getTwoFragmentTargets();

  if (!t) {
    twoResultNg.textContent = "—";
    twoResultDetail.textContent = "";
    updateTwoMasterMix();
    return;
  }

  if (twoKnow === "vector") {
    twoResultNg.textContent = `${fmt(t.insertNg)} ng`;
    twoResultDetail.textContent =
      `${fmt(t.vectorPmol, 4)} pmol vector → ${fmt(t.insertPmol, 4)} pmol insert needed (${fmt(t.ratio)}:1)`;
  } else {
    twoResultNg.textContent = `${fmt(t.vectorNg)} ng`;
    twoResultDetail.textContent =
      `${fmt(t.insertPmol, 4)} pmol insert → ${fmt(t.vectorPmol, 4)} pmol vector needed (ratio ${fmt(t.ratio)}:1)`;
  }

  updateTwoMasterMix();
}

[twoVectorLen, twoVectorUnit, twoInsertLen, twoInsertUnit, twoRatio, twoKnownNg].forEach((el) =>
  el.addEventListener("input", calcTwoFragment)
);

// ---- 2-fragment master mix ----
const twoMmVectorConc = document.getElementById("two-mm-vector-conc");
const twoMmInsertConc = document.getElementById("two-mm-insert-conc");
const twoMmVectorVol = document.getElementById("two-mm-vector-vol");
const twoMmInsertVol = document.getElementById("two-mm-insert-vol");
const twoMmWater = document.getElementById("two-mm-water");
const twoMmMm = document.getElementById("two-mm-mm");
const twoMmTotal = document.getElementById("two-mm-total");
const twoMmStatus = document.getElementById("two-mm-status");
const twoMmDnaVol = document.getElementById("two-mm-dnavol");
const twoMmIdeal = document.getElementById("two-mm-ideal");
const twoMmFallback = document.getElementById("two-mm-fallback");

function updateTwoMasterMix() {
  const t = getTwoFragmentTargets();
  const vectorConc = parseFloat(twoMmVectorConc.value);
  const insertConc = parseFloat(twoMmInsertConc.value);
  const dnaVol = parseFloat(twoMmDnaVol.value);
  const idealMin = parseFloat(twoMmIdeal.value);
  const fallbackMin = parseFloat(twoMmFallback.value);

  const ready = t && vectorConc > 0 && insertConc > 0 && dnaVol > 0 && idealMin > 0 && fallbackMin > 0 && fallbackMin <= idealMin;

  if (!ready) {
    twoMmVectorVol.textContent = "—";
    twoMmInsertVol.textContent = "—";
    twoMmWater.textContent = "—";
    twoMmMm.textContent = "—";
    twoMmTotal.textContent = "—";
    twoMmStatus.textContent = "";
    return;
  }

  const result = solveMasterMix(
    [
      { label: "vector", targetNg: t.vectorNg, stockConc: vectorConc },
      { label: "insert", targetNg: t.insertNg, stockConc: insertConc },
    ],
    dnaVol,
    idealMin,
    fallbackMin
  );

  renderMmVolCell(twoMmVectorVol, result.rows[0]);
  renderMmVolCell(twoMmInsertVol, result.rows[1]);
  twoMmWater.textContent = `${fmt(result.water)} µL`;
  twoMmMm.textContent = `${fmt(result.mmVol)} µL`;
  twoMmTotal.textContent = `${fmt(result.totalReaction)} µL`;
  twoMmStatus.textContent = masterMixStatusText(result, idealMin, fallbackMin);
}

[twoMmVectorConc, twoMmInsertConc, twoMmDnaVol, twoMmIdeal, twoMmFallback].forEach((el) =>
  el.addEventListener("input", updateTwoMasterMix)
);

// ================= MULTI-FRAGMENT MODE =================
const multiPmol = document.getElementById("multi-pmol");
const fragmentRows = document.getElementById("fragment-rows");
const addFragmentBtn = document.getElementById("add-fragment");
const multiResultTotal = document.getElementById("multi-result-total");

let fragmentCount = 0;

function addFragmentRow(name = "", len = "", unit = "bp") {
  fragmentCount += 1;
  const row = document.createElement("tr");
  row.innerHTML = `
    <td><input type="text" class="frag-name" value="${name || `Fragment ${fragmentCount}`}"></td>
    <td>
      <div class="length-input">
        <input type="number" class="frag-len" min="0" step="any" value="${len}">
        <select class="frag-unit">
          <option value="bp" ${unit === "bp" ? "selected" : ""}>bp</option>
          <option value="kb" ${unit === "kb" ? "selected" : ""}>kb</option>
        </select>
      </div>
    </td>
    <td class="fragment-ng">—</td>
    <td><input type="number" class="frag-conc" min="0" step="any" placeholder="ng/µL"></td>
    <td class="mm-vol frag-pipette">—</td>
    <td><button class="remove-row-btn" title="Remove fragment">✕</button></td>
  `;
  fragmentRows.appendChild(row);

  row.querySelector(".frag-len").addEventListener("input", calcMultiFragment);
  row.querySelector(".frag-unit").addEventListener("change", calcMultiFragment);
  row.querySelector(".frag-conc").addEventListener("input", calcMultiFragment);
  row.querySelector(".remove-row-btn").addEventListener("click", () => {
    row.remove();
    calcMultiFragment();
  });

  calcMultiFragment();
}

function calcMultiFragment() {
  const pmol = parseFloat(multiPmol.value);
  let total = 0;
  let anyValid = false;

  fragmentRows.querySelectorAll("tr").forEach((row) => {
    const bp = bpFromLength(row.querySelector(".frag-len").value, row.querySelector(".frag-unit").value);
    const cell = row.querySelector(".fragment-ng");
    if (!pmol || !bp || bp <= 0) {
      cell.textContent = "—";
      return;
    }
    const ng = ngFromPmol(pmol, bp);
    cell.textContent = `${fmt(ng)} ng`;
    total += ng;
    anyValid = true;
  });

  multiResultTotal.textContent = anyValid ? `${fmt(total)} ng` : "—";
  updateMultiMasterMix();
}

addFragmentBtn.addEventListener("click", () => addFragmentRow());
multiPmol.addEventListener("input", calcMultiFragment);

// ---- multi-fragment master mix ----
const multiMmWater = document.getElementById("multi-mm-water");
const multiMmMm = document.getElementById("multi-mm-mm");
const multiMmTotal = document.getElementById("multi-mm-total");
const multiMmStatus = document.getElementById("multi-mm-status");
const multiMmDnaVol = document.getElementById("multi-mm-dnavol");
const multiMmIdeal = document.getElementById("multi-mm-ideal");
const multiMmFallback = document.getElementById("multi-mm-fallback");

function updateMultiMasterMix() {
  const pmol = parseFloat(multiPmol.value);
  const dnaVol = parseFloat(multiMmDnaVol.value);
  const idealMin = parseFloat(multiMmIdeal.value);
  const fallbackMin = parseFloat(multiMmFallback.value);

  const rowEls = Array.from(fragmentRows.querySelectorAll("tr")).filter((row) => {
    const bp = bpFromLength(row.querySelector(".frag-len").value, row.querySelector(".frag-unit").value);
    return !!bp && bp > 0;
  });

  const settingsValid = pmol > 0 && dnaVol > 0 && idealMin > 0 && fallbackMin > 0 && fallbackMin <= idealMin;

  const fragments = rowEls.map((row) => {
    const bp = bpFromLength(row.querySelector(".frag-len").value, row.querySelector(".frag-unit").value);
    const stockConc = parseFloat(row.querySelector(".frag-conc").value);
    return { row, targetNg: ngFromPmol(pmol, bp), stockConc };
  });

  const allHaveConc = fragments.length > 0 && fragments.every((f) => f.stockConc > 0);

  if (!settingsValid || !allHaveConc) {
    fragments.forEach((f) => {
      f.row.querySelector(".frag-pipette").textContent = "—";
    });
    multiMmWater.textContent = "—";
    multiMmMm.textContent = "—";
    multiMmTotal.textContent = "—";
    multiMmStatus.textContent = fragments.length > 0 ? "Enter a stock concentration for every fragment." : "";
    return;
  }

  const result = solveMasterMix(fragments, dnaVol, idealMin, fallbackMin);

  result.rows.forEach((r) => renderMmVolCell(r.row.querySelector(".frag-pipette"), r));

  multiMmWater.textContent = `${fmt(result.water)} µL`;
  multiMmMm.textContent = `${fmt(result.mmVol)} µL`;
  multiMmTotal.textContent = `${fmt(result.totalReaction)} µL`;
  multiMmStatus.textContent = masterMixStatusText(result, idealMin, fallbackMin);
}

[multiMmDnaVol, multiMmIdeal, multiMmFallback].forEach((el) =>
  el.addEventListener("input", updateMultiMasterMix)
);

// seed with a typical 3-fragment example (vector + 2 inserts)
addFragmentRow("Vector", 3000);
addFragmentRow("Insert 1", 1000);
addFragmentRow("Insert 2", 800);

// initial calc
calcTwoFragment();

// ================= CRISPR ELECTROPORATION MIX =================
const crGuideLen = document.getElementById("cr-guide-len");
const crGuideUnit = document.getElementById("cr-guide-unit");
const crGuideConc = document.getElementById("cr-guide-conc");
const crGuideVol = document.getElementById("cr-guide-vol");
const crInsertLen = document.getElementById("cr-insert-len");
const crInsertUnit = document.getElementById("cr-insert-unit");
const crInsertConc = document.getElementById("cr-insert-conc");
const crInsertVol = document.getElementById("cr-insert-vol");
const crGuideParts = document.getElementById("cr-guide-parts");
const crInsertParts = document.getElementById("cr-insert-parts");
const crGuideResult = document.getElementById("cr-guide-result");
const crGuideDetail = document.getElementById("cr-guide-detail");
const crInsertResult = document.getElementById("cr-insert-result");
const crInsertDetail = document.getElementById("cr-insert-detail");
const crLimiting = document.getElementById("cr-limiting");
const crTotalDetail = document.getElementById("cr-total-detail");

const CR_LOW_VOL_WARNING = 1; // µL — below this, flag that dilution may help

function calcCrispr() {
  const guideBp = bpFromLength(crGuideLen.value, crGuideUnit.value);
  const insertBp = bpFromLength(crInsertLen.value, crInsertUnit.value);
  const guideConc = parseFloat(crGuideConc.value);
  const insertConc = parseFloat(crInsertConc.value);
  const guideMaxVol = parseFloat(crGuideVol.value);
  const insertMaxVol = parseFloat(crInsertVol.value);
  const guideParts = parseFloat(crGuideParts.value);
  const insertParts = parseFloat(crInsertParts.value);

  const valid =
    guideBp && insertBp && guideConc > 0 && insertConc > 0 &&
    guideMaxVol > 0 && insertMaxVol > 0 && guideParts > 0 && insertParts > 0;

  if (!valid) {
    crGuideResult.textContent = "—";
    crGuideDetail.textContent = "";
    crInsertResult.textContent = "—";
    crInsertDetail.textContent = "";
    crLimiting.textContent = "—";
    crTotalDetail.textContent = "";
    return;
  }

  // limiting-reagent stoichiometry: whichever prep yields fewer pmol per ratio part caps the mix
  const guideMaxPmol = pmolFromNg(guideConc * guideMaxVol, guideBp);
  const insertMaxPmol = pmolFromNg(insertConc * insertMaxVol, insertBp);

  const guideUnitPmol = guideMaxPmol / guideParts;
  const insertUnitPmol = insertMaxPmol / insertParts;
  const unitPmol = Math.min(guideUnitPmol, insertUnitPmol);
  const limiting = guideUnitPmol <= insertUnitPmol ? "guide" : "insert";

  const guidePmolUsed = unitPmol * guideParts;
  const insertPmolUsed = unitPmol * insertParts;

  const guideNgUsed = ngFromPmol(guidePmolUsed, guideBp);
  const insertNgUsed = ngFromPmol(insertPmolUsed, insertBp);

  const guideVolUsed = guideNgUsed / guideConc;
  const insertVolUsed = insertNgUsed / insertConc;

  const lowVolNote = (vol) =>
    vol < CR_LOW_VOL_WARNING ? '<span class="mm-note">⚠ below 1 µL — consider diluting stock</span>' : "";

  crGuideResult.innerHTML = `${fmt(guideVolUsed)} µL${lowVolNote(guideVolUsed)}`;
  crGuideDetail.textContent = `${fmt(guideNgUsed)} ng · ${fmt(guidePmolUsed, 4)} pmol`;

  crInsertResult.innerHTML = `${fmt(insertVolUsed)} µL${lowVolNote(insertVolUsed)}`;
  crInsertDetail.textContent = `${fmt(insertNgUsed)} ng · ${fmt(insertPmolUsed, 4)} pmol`;

  crLimiting.textContent = limiting === "guide" ? "Guide" : "Insert";
  crTotalDetail.textContent =
    `${limiting === "guide" ? "Guide" : "Insert"} uses its full available volume ` +
    `(${fmt(limiting === "guide" ? guideMaxVol : insertMaxVol)} µL). ` +
    `Combined DNA volume: ${fmt(guideVolUsed + insertVolUsed)} µL.`;
}

[
  crGuideLen, crGuideUnit, crGuideConc, crGuideVol,
  crInsertLen, crInsertUnit, crInsertConc, crInsertVol,
  crGuideParts, crInsertParts,
].forEach((el) => el.addEventListener("input", calcCrispr));

calcCrispr();

// ================= WESTERN BLOT SAMPLE BUFFER =================
const wbBufferFactor = document.getElementById("wb-buffer-factor");
const wbDttFactor = document.getElementById("wb-dtt-factor");
const wbKnowSwitch = document.getElementById("wb-know-switch");
const wbKnownTitle = document.getElementById("wb-known-title");
const wbKnownLabel = document.getElementById("wb-known-label");
const wbKnownVol = document.getElementById("wb-known-vol");
const wbPresets = document.getElementById("wb-presets");
const wbSampleVol = document.getElementById("wb-sample-vol");
const wbDttVol = document.getElementById("wb-dtt-vol");
const wbBufferVol = document.getElementById("wb-buffer-vol");
const wbTotalVol = document.getElementById("wb-total-vol");
const wbStatus = document.getElementById("wb-status");
const wbFormula = document.getElementById("wb-formula");

let wbKnow = "total"; // "total" = have total volume, solve components. "sample" = have sample volume.

wbKnowSwitch.querySelectorAll(".mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    wbKnowSwitch.querySelectorAll(".mode-btn").forEach((b) => b.classList.remove("active"));
    btn.classList.add("active");
    wbKnow = btn.dataset.know;

    if (wbKnow === "total") {
      wbKnownTitle.textContent = "Total volume";
      wbPresets.style.display = "";
      wbFormula.textContent = "buffer = total ÷ buffer factor; DTT = total ÷ DTT factor; sample = total − buffer − DTT";
    } else {
      wbKnownTitle.textContent = "Sample volume";
      wbPresets.style.display = "none";
      wbFormula.textContent = "total = sample ÷ (1 − 1/buffer factor − 1/DTT factor); buffer = total ÷ buffer factor; DTT = total ÷ DTT factor";
    }
    calcWestern();
  });
});

wbPresets.querySelectorAll(".preset-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    wbKnownVol.value = btn.dataset.preset;
    calcWestern();
  });
});

function calcWestern() {
  const bufferFactor = parseFloat(wbBufferFactor.value);
  const dttFactor = parseFloat(wbDttFactor.value);
  const knownVol = parseFloat(wbKnownVol.value);

  const validFactors = bufferFactor > 1 && dttFactor > 1;
  const sampleFraction = validFactors ? 1 - 1 / bufferFactor - 1 / dttFactor : NaN;

  if (!validFactors || sampleFraction <= 0 || !knownVol || knownVol <= 0) {
    [wbSampleVol, wbDttVol, wbBufferVol, wbTotalVol].forEach((el) => (el.textContent = "—"));
    wbStatus.textContent =
      validFactors && sampleFraction <= 0
        ? "Buffer + DTT factors leave no volume for sample — use larger dilution factors."
        : "";
    return;
  }

  const total = wbKnow === "total" ? knownVol : knownVol / sampleFraction;
  const buffer = total / bufferFactor;
  const dtt = total / dttFactor;
  const sample = total - buffer - dtt;

  wbSampleVol.textContent = `${fmt(sample)} µL`;
  wbDttVol.textContent = `${fmt(dtt)} µL`;
  wbBufferVol.textContent = `${fmt(buffer)} µL`;
  wbTotalVol.textContent = `${fmt(total)} µL`;
  wbStatus.textContent = "";
}

[wbBufferFactor, wbDttFactor, wbKnownVol].forEach((el) => el.addEventListener("input", calcWestern));

calcWestern();

// ================= PRIMER Tm / Ta =================
// Replicates the NEB Tm Calculator's method: SantaLucia (1998) unified nearest-neighbor
// thermodynamics, Owczarzy (2004) / Schildkraut (1965) salt correction, and each enzyme's
// published empirical Ta adjustment. Parameters below were reverse-engineered from NEB's own
// public client-side calculator (tmcalculator.neb.com) for the enzymes listed.
const NN_SANTALUCIA = {
  aa: { dh: -7.9, ds: -22.2 }, tt: { dh: -7.9, ds: -22.2 },
  at: { dh: -7.2, ds: -20.4 }, ta: { dh: -7.2, ds: -21.3 },
  ca: { dh: -8.5, ds: -22.7 }, tg: { dh: -8.5, ds: -22.7 },
  gt: { dh: -8.4, ds: -22.4 }, ac: { dh: -8.4, ds: -22.4 },
  ct: { dh: -7.8, ds: -21 }, ag: { dh: -7.8, ds: -21 },
  ga: { dh: -8.2, ds: -22.2 }, tc: { dh: -8.2, ds: -22.2 },
  cg: { dh: -10.6, ds: -27.2 }, gc: { dh: -9.8, ds: -24.4 },
  gg: { dh: -8, ds: -19.9 }, cc: { dh: -8, ds: -19.9 },
};

// method 4 = Owczarzy salt correction (most enzymes); method 5 = Schildkraut + Ct/4 (Phusion)
const POLY_DATA = {
  q5: { label: "Q5 High-Fidelity", saltMM: 150, ctNM: 500, method: 4, taRule: "q5plus1" },
  q5hs: { label: "Q5 Hot Start", saltMM: 150, ctNM: 500, method: 4, taRule: "q5plus1" },
  q5bd: { label: "Q5 Blood Direct", saltMM: 150, ctNM: 500, method: 4, taRule: "q5plus1" },
  q5u: { label: "Q5U Hot Start", saltMM: 170, ctNM: 500, method: 4, taRule: "q5uplus2" },
  phusion: { label: "Phusion", saltMM: 222, ctNM: 500, method: 5, taRule: "phusion" },
  phusionflex: { label: "Phusion Hot Start Flex", saltMM: 222, ctNM: 500, method: 5, taRule: "phusion" },
  hstaq: { label: "Hot Start Taq", saltMM: 55, ctNM: 200, method: 4, taRule: "taqlike68" },
  taq: { label: "Taq DNA Polymerase", saltMM: 55, ctNM: 200, method: 4, taRule: "taqlike68" },
  onetaq: { label: "OneTaq (Standard Buffer)", saltMM: 54, ctNM: 200, method: 4, taRule: "taqlike68" },
  onetaqhs: { label: "OneTaq Hot Start (Standard Buffer)", saltMM: 54, ctNM: 200, method: 4, taRule: "taqlike68" },
  hemoklentaq: { label: "Hemo KlenTaq", saltMM: 70, ctNM: 200, method: 4, taRule: "taqlike68" },
  longamp: { label: "LongAmp Taq", saltMM: 100, ctNM: 400, method: 4, taRule: "taqlike65" },
  longamphs: { label: "LongAmp Hot Start Taq", saltMM: 100, ctNM: 400, method: 4, taRule: "taqlike65" },
  custom: { label: "Custom (enter salt below)", saltMM: 50, ctNM: 250, method: 4, taRule: null },
};

function prRevcomp(seq) {
  const comp = { a: "t", t: "a", g: "c", c: "g" };
  return seq.split("").reverse().map((ch) => comp[ch]).join("");
}

function prNormalize(rawSeq) {
  const seq = (rawSeq || "").toLowerCase().replace(/\s/g, "").replace(/u/g, "t");
  return /^[acgt]+$/.test(seq) ? seq : null;
}

// nearest-neighbor enthalpy/entropy sum for a short duplex formed by `seq` binding its own
// perfect complement — used both for full-primer Tm and for scoring individual dimer runs
function nnThermo(seq) {
  const len = seq.length;
  let dh = 0, ds = 0; // NN sums: dh in kcal/mol, ds in cal/mol/K
  for (let i = 0; i < len - 1; i++) {
    const nn = NN_SANTALUCIA[seq.slice(i, i + 2)];
    dh += nn.dh;
    ds += nn.ds;
  }

  let initDs = 0, initDhCal = 0;
  for (const ch of [seq[0], seq[len - 1]]) {
    if (ch === "a" || ch === "t") { initDs += 4.1; initDhCal += 2300; }
    else { initDs += -2.8; initDhCal += 100; }
  }
  const sym = len > 1 && seq === prRevcomp(seq);
  const symDs = sym ? -1.4 : 0;

  return { dhCal: dh * 1000 + initDhCal, ds: ds + initDs + symDs, sym };
}

// ΔG (kcal/mol) of `seq` binding its own perfect complement at a fixed reference temperature.
// Uses raw (1 M NaCl reference) SantaLucia parameters, unadjusted for reaction salt — the
// standard way oligo ΔG is reported, since IDT's actual salt handling for this value isn't public.
function duplexDeltaG(seq, tempC = 25) {
  if (!seq || seq.length < 2) return null;
  const { dhCal, ds } = nnThermo(seq);
  const tK = tempC + 273.15;
  return (dhCal - tK * ds) / 1000;
}

// core NN Tm calculation; seq is raw user input, saltMM/ctNM/method describe the reaction buffer
function calcPrimerTm(rawSeq, saltMM, ctNM, method) {
  const seq = prNormalize(rawSeq);
  if (!seq || seq.length < 2) return null;

  const len = seq.length;
  let gcCount = 0;
  for (const ch of seq) if (ch === "g" || ch === "c") gcCount++;
  const fgc = gcCount / len;

  const { dhCal: totalDhCal, ds: totalDs } = nnThermo(seq);

  const R = 1.987;
  const ctM = (method === 5 ? ctNM / 4 : ctNM) * 1e-9;
  const saltM = saltMM / 1000;

  const tmK = totalDhCal / (totalDs + R * Math.log(ctM));

  let correctedK;
  if (method === 5) {
    correctedK = tmK + 16.6 * (Math.log(saltM) / Math.LN10);
  } else {
    const scOw = 1e-5 * (4.29 * fgc - 3.95) * Math.log(saltM) + 9.4e-6 * Math.pow(Math.log(saltM), 2);
    correctedK = 1 / (1 / tmK + scOw);
  }

  // NEB rounds Tm to the nearest whole degree before using it in the Ta formula (and for display)
  return { len, gc: Math.round(fgc * 1000) / 10, tm: Math.round(correctedK - 273.15) };
}

// Scans every possible antiparallel offset between seqA (5'->3') and seqB (5'->3', reversed to
// lie 3'->5' underneath it) for the longest run of Watson-Crick complementary bases — the same
// "longest stretch of complementary bases" approach IDT's OligoAnalyzer describes for its
// self-dimer/hetero-dimer check. Returns up to `maxResults` candidate dimers, strongest ΔG first.
function findDimerRegions(rawA, rawB, minRun = 2, maxResults = 3) {
  const seqA = prNormalize(rawA);
  const seqB = prNormalize(rawB);
  if (!seqA || !seqB) return null;

  const comp = { a: "t", t: "a", g: "c", c: "g" };
  const revB = seqB.split("").reverse().join("");
  const nA = seqA.length, nB = revB.length;

  const candidates = [];
  for (let shift = -(nB - 1); shift <= nA - 1; shift++) {
    let curStart = -1, curLen = 0, bestStart = -1, bestLen = 0;
    const lo = Math.max(0, shift), hi = Math.min(nA - 1, nB - 1 + shift);
    for (let i = lo; i <= hi + 1; i++) {
      const match = i <= hi && comp[seqA[i]] === revB[i - shift];
      if (match) {
        if (curLen === 0) curStart = i;
        curLen++;
      } else {
        if (curLen > bestLen) { bestLen = curLen; bestStart = curStart; }
        curLen = 0;
      }
    }
    if (bestLen >= minRun) {
      const topRun = seqA.slice(bestStart, bestStart + bestLen);
      candidates.push({ shift, start: bestStart, len: bestLen, dG: duplexDeltaG(topRun) });
    }
  }

  candidates.sort((a, b) => a.dG - b.dG);

  // collapse near-duplicate registers (adjacent shifts often rediscover the same core run)
  const deduped = [];
  for (const c of candidates) {
    if (deduped.some((d) => Math.abs(d.shift - c.shift) <= 1 && Math.abs(d.dG - c.dG) < 0.01)) continue;
    deduped.push(c);
    if (deduped.length >= maxResults) break;
  }

  return deduped.map((c) => ({ ...c, align: renderDimerAlignment(seqA, seqB, c.shift, c.start, c.len) }));
}

// builds a 3-line monospace alignment (top strand / pairing marks / bottom strand) for one register
function renderDimerAlignment(seqA, seqB, shift, start, len) {
  const revB = seqB.split("").reverse().join("");
  const padTop = shift < 0 ? -shift : 0;
  const padBottom = shift > 0 ? shift : 0;
  const top = " ".repeat(padTop) + seqA.toUpperCase();
  const bottom = " ".repeat(padBottom) + revB.toUpperCase();
  const width = Math.max(top.length, bottom.length);
  const pairLine = Array(width).fill(" ");
  for (let i = start; i < start + len; i++) pairLine[i + padTop] = "|";
  return {
    top: top.padEnd(width),
    pair: pairLine.join(""),
    bottom: bottom.padEnd(width),
  };
}

function computeTa(rule, minTm, minLen) {
  let a = minTm;
  switch (rule) {
    case "q5plus1":
      if (minLen > 7) a = minTm + 1;
      if (a > 72) a = 72;
      break;
    case "q5uplus2":
      if (minLen > 7) a = minTm + 2;
      if (a > 72) a = 72;
      break;
    case "phusion":
      a = 0.93 * minTm + 7.5;
      if (a > 72) a = 72;
      break;
    case "taqlike68":
      if (minLen > 7) a = minTm - 5;
      if (a > 68) a = 68;
      break;
    case "taqlike65":
      if (minLen > 7) a = minTm - 5;
      if (a > 65) a = 65;
      break;
    default:
      return null;
  }
  return Math.round(a * 10) / 10;
}

// clipboard helper with a fallback for contexts where the Clipboard API is unavailable (e.g. file://)
function copyText(text) {
  if (navigator.clipboard && window.isSecureContext) {
    return navigator.clipboard.writeText(text);
  }
  const ta = document.createElement("textarea");
  ta.value = text;
  ta.style.position = "fixed";
  ta.style.opacity = "0";
  document.body.appendChild(ta);
  ta.focus();
  ta.select();
  try {
    document.execCommand("copy");
  } finally {
    document.body.removeChild(ta);
  }
  return Promise.resolve();
}

const prEnzyme = document.getElementById("pr-enzyme");
const prConc = document.getElementById("pr-conc");
const prSaltWrap = document.getElementById("pr-salt-wrap");
const prSalt = document.getElementById("pr-salt");
const prSeq1 = document.getElementById("pr-seq1");
const prSeq2 = document.getElementById("pr-seq2");
const prLen1 = document.getElementById("pr-len1");
const prGc1 = document.getElementById("pr-gc1");
const prTm1 = document.getElementById("pr-tm1");
const prLen2 = document.getElementById("pr-len2");
const prGc2 = document.getElementById("pr-gc2");
const prTm2 = document.getElementById("pr-tm2");
const prTa = document.getElementById("pr-ta");
const prStatus = document.getElementById("pr-status");
const prCopy1 = document.getElementById("pr-copy1");
const prCopy2 = document.getElementById("pr-copy2");

Object.keys(POLY_DATA).forEach((key) => {
  const opt = document.createElement("option");
  opt.value = key;
  opt.textContent = POLY_DATA[key].label;
  if (key === "q5") opt.selected = true;
  prEnzyme.appendChild(opt);
});

prEnzyme.addEventListener("change", () => {
  const data = POLY_DATA[prEnzyme.value];
  prConc.value = data.ctNM;
  prSaltWrap.style.display = prEnzyme.value === "custom" ? "" : "none";
  calcPrimers();
});

function calcPrimers() {
  const data = POLY_DATA[prEnzyme.value];
  const ctNM = parseFloat(prConc.value);
  const saltMM = prEnzyme.value === "custom" ? parseFloat(prSalt.value) : data.saltMM;

  if (!ctNM || ctNM <= 0 || !saltMM || saltMM <= 0) {
    [prLen1, prGc1, prTm1, prLen2, prGc2, prTm2, prTa].forEach((el) => (el.textContent = "—"));
    prStatus.textContent = "";
    calcDimers();
    return;
  }

  const seq1raw = prSeq1.value.trim();
  const seq2raw = prSeq2.value.trim();
  const r1 = seq1raw ? calcPrimerTm(seq1raw, saltMM, ctNM, data.method) : null;
  const r2 = seq2raw ? calcPrimerTm(seq2raw, saltMM, ctNM, data.method) : null;

  prLen1.textContent = r1 ? r1.len : "—";
  prGc1.textContent = r1 ? `${fmt(r1.gc)}%` : "—";
  prTm1.textContent = r1 ? `${fmt(r1.tm)}°C` : "—";

  prLen2.textContent = r2 ? r2.len : "—";
  prGc2.textContent = r2 ? `${fmt(r2.gc)}%` : "—";
  prTm2.textContent = r2 ? `${fmt(r2.tm)}°C` : "—";

  let status = "";
  if (seq1raw && !r1) status = "Primer 1: only A/C/G/T/U characters are supported.";
  else if (seq2raw && !r2) status = "Primer 2: only A/C/G/T/U characters are supported.";

  if (r1 && r2 && data.taRule) {
    const minTm = Math.min(r1.tm, r2.tm);
    const minLen = Math.min(r1.len, r2.len);
    prTa.textContent = `${fmt(computeTa(data.taRule, minTm, minLen))}°C`;
  } else {
    prTa.textContent = "—";
    if (r1 && r2 && !data.taRule) status = "Ta isn't available for Custom — pick a listed enzyme for an annealing temp.";
  }

  prStatus.textContent = status;

  calcDimers();
}

const prDimerSelf1 = document.getElementById("pr-dimer-self1");
const prDimerSelf2 = document.getElementById("pr-dimer-self2");
const prDimerHetero = document.getElementById("pr-dimer-hetero");

function renderDimerBlock(container, seqA, seqB, maxDeltaGSeq) {
  const seq = prNormalize(seqA);
  const seqOther = prNormalize(seqB);
  if (!seq || !seqOther) {
    container.innerHTML = `<div class="dimer-none">Enter a valid sequence to check.</div>`;
    return;
  }

  const regions = findDimerRegions(seqA, seqB);
  let html = "";

  if (maxDeltaGSeq) {
    const maxDg = duplexDeltaG(prNormalize(maxDeltaGSeq));
    html += `<div class="dimer-none">Maximum ΔG (binding its perfect complement): ${fmt(maxDg)} kcal/mol</div>`;
  }

  if (!regions || regions.length === 0) {
    html += `<div class="dimer-none">No self-complementary run of 2+ bp found.</div>`;
  } else {
    regions.forEach((r) => {
      html += `
        <div class="dimer-result">
          <div class="dimer-result-header">
            <span class="dimer-dg">ΔG: ${fmt(r.dG)} kcal/mol</span>
            <span class="table-note">${r.len} bp</span>
          </div>
          <pre class="dimer-align">5' ${r.align.top} 3'
   ${r.align.pair}
3' ${r.align.bottom} 5'</pre>
        </div>`;
    });
  }

  container.innerHTML = html;
}

function calcDimers() {
  const seq1raw = prSeq1.value.trim();
  const seq2raw = prSeq2.value.trim();

  renderDimerBlock(prDimerSelf1, seq1raw, seq1raw, seq1raw);
  renderDimerBlock(prDimerSelf2, seq2raw, seq2raw, seq2raw);

  if (seq1raw && seq2raw) {
    renderDimerBlock(prDimerHetero, seq1raw, seq2raw, null);
  } else {
    prDimerHetero.innerHTML = `<div class="dimer-none">Enter both primers to check.</div>`;
  }
}

function wireCopyButton(button, seqInput, url, copiedText) {
  button.addEventListener("click", () => {
    const seq = seqInput.value.trim();
    if (!seq) return;
    copyText(seq).then(() => {
      const original = button.textContent;
      button.textContent = copiedText;
      window.open(url, "_blank", "noopener");
      setTimeout(() => (button.textContent = original), 3000);
    });
  });
}

const MFOLD_URL = "https://www.unafold.org/mfold/applications/dna-folding-form.php";
const OLIGOANALYZER_URL = "https://www.idtdna.com/pages/tools/oligoanalyzer";

wireCopyButton(document.getElementById("pr-copy1-mfold"), prSeq1, MFOLD_URL, "Copied! Paste into the Sequence box →");
wireCopyButton(document.getElementById("pr-copy2-mfold"), prSeq2, MFOLD_URL, "Copied! Paste into the Sequence box →");
wireCopyButton(prCopy1, prSeq1, OLIGOANALYZER_URL, "Copied! Sign in, then paste it in →");
wireCopyButton(prCopy2, prSeq2, OLIGOANALYZER_URL, "Copied! Sign in, then paste it in →");

[prConc, prSalt, prSeq1, prSeq2].forEach((el) => el.addEventListener("input", calcPrimers));

calcPrimers();

// ================= EXPORT (.txt) =================
// Each builder recomputes fresh from the raw input state (reusing the same functions the
// on-screen calculators use) rather than scraping rendered DOM/HTML, so the exported text is
// always clean plain text even where the on-screen cell mixes a value with an inline note.
function downloadText(filename, text) {
  const blob = new Blob([text], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function exportFilename(tool) {
  const stamp = new Date().toISOString().slice(0, 10);
  return `${tool}-${stamp}.txt`;
}

function buildHifiExport() {
  const lines = ["HiFi (Gibson) Assembly DNA Ratio Calculator", `Generated: ${new Date().toLocaleString()}`, ""];
  const twoActive = document.getElementById("mode-two").classList.contains("active");

  if (twoActive) {
    lines.push("Mode: 2-Fragment (Vector + Insert)", "");
    lines.push(`Vector length: ${twoVectorLen.value} ${twoVectorUnit.value}`);
    lines.push(`Insert length: ${twoInsertLen.value} ${twoInsertUnit.value}`);
    lines.push(`Ratio (insert:vector): ${twoRatio.value}`);
    lines.push(`Known: ${twoKnownTitle.textContent} = ${twoKnownNg.value} ng`, "");

    const t = getTwoFragmentTargets();
    if (t) {
      lines.push(`${twoResultLabel.textContent}: ${twoResultNg.textContent}`);
      lines.push(`  ${twoResultDetail.textContent}`);
    } else {
      lines.push("Result: (enter valid values)");
    }

    lines.push("", "Master Mix Setup");
    const vectorConc = parseFloat(twoMmVectorConc.value);
    const insertConc = parseFloat(twoMmInsertConc.value);
    const dnaVol = parseFloat(twoMmDnaVol.value);
    const idealMin = parseFloat(twoMmIdeal.value);
    const fallbackMin = parseFloat(twoMmFallback.value);
    const ready = t && vectorConc > 0 && insertConc > 0 && dnaVol > 0 && idealMin > 0 && fallbackMin > 0 && fallbackMin <= idealMin;

    if (ready) {
      const result = solveMasterMix(
        [{ targetNg: t.vectorNg, stockConc: vectorConc }, { targetNg: t.insertNg, stockConc: insertConc }],
        dnaVol, idealMin, fallbackMin
      );
      const rowLine = (label, row) => {
        let s = `  ${label}: ${fmt(row.vol)} µL`;
        if (row.diluted) s += ` (dilute to ${fmt(row.dilutedConc)} ng/µL, ratio 1:${fmt(row.factor)})`;
        return s;
      };
      lines.push(rowLine("Vector", result.rows[0]));
      lines.push(rowLine("Insert", result.rows[1]));
      lines.push(`  ddH2O: ${fmt(result.water)} µL`);
      lines.push(`  2x HiFi MM: ${fmt(result.mmVol)} µL`);
      lines.push(`  Total reaction: ${fmt(result.totalReaction)} µL`);
      lines.push(`  ${masterMixStatusText(result, idealMin, fallbackMin)}`);
    } else {
      lines.push("  (enter stock concentrations to compute)");
    }
  } else {
    lines.push("Mode: Multi-Fragment (Equimolar)", "");
    const pmol = parseFloat(multiPmol.value);
    lines.push(`Target per fragment: ${multiPmol.value} pmol`, "");

    const rowEls = Array.from(fragmentRows.querySelectorAll("tr"));
    const fragments = rowEls.map((row) => {
      const name = row.querySelector(".frag-name").value || "Fragment";
      const bp = bpFromLength(row.querySelector(".frag-len").value, row.querySelector(".frag-unit").value);
      const stockConc = parseFloat(row.querySelector(".frag-conc").value);
      return { name, bp, stockConc, targetNg: bp ? ngFromPmol(pmol, bp) : NaN };
    });

    let totalNg = 0;
    const validFrags = fragments.filter((f) => f.bp);
    validFrags.forEach((f) => { totalNg += f.targetNg; });

    const dnaVol = parseFloat(multiMmDnaVol.value);
    const idealMin = parseFloat(multiMmIdeal.value);
    const fallbackMin = parseFloat(multiMmFallback.value);
    const allHaveConc = validFrags.length > 0 && validFrags.every((f) => f.stockConc > 0);
    const settingsValid = pmol > 0 && dnaVol > 0 && idealMin > 0 && fallbackMin > 0 && fallbackMin <= idealMin;

    let mmResult = null;
    if (settingsValid && allHaveConc) {
      mmResult = solveMasterMix(validFrags.map((f) => ({ targetNg: f.targetNg, stockConc: f.stockConc })), dnaVol, idealMin, fallbackMin);
    }

    fragments.forEach((f) => {
      if (!f.bp) { lines.push(`${f.name}: (enter length)`); return; }
      let line = `${f.name}: ${fmt(f.bp)} bp — target ${fmt(f.targetNg)} ng`;
      if (f.stockConc > 0) line += `, stock ${fmt(f.stockConc)} ng/µL`;
      lines.push(line);
      if (mmResult) {
        const idx = validFrags.indexOf(f);
        const row = mmResult.rows[idx];
        let pipetteLine = `  Pipette: ${fmt(row.vol)} µL`;
        if (row.diluted) pipetteLine += ` (dilute to ${fmt(row.dilutedConc)} ng/µL, ratio 1:${fmt(row.factor)})`;
        lines.push(pipetteLine);
      }
    });

    lines.push("", `Total DNA (by mass): ${fmt(totalNg)} ng`, "", "Master Mix Setup");
    if (mmResult) {
      lines.push(`  ddH2O: ${fmt(mmResult.water)} µL`);
      lines.push(`  2x HiFi MM: ${fmt(mmResult.mmVol)} µL`);
      lines.push(`  Total reaction: ${fmt(mmResult.totalReaction)} µL`);
      lines.push(`  ${masterMixStatusText(mmResult, idealMin, fallbackMin)}`);
    } else {
      lines.push("  (enter stock concentrations for every fragment to compute)");
    }
  }

  return lines.join("\n");
}

function buildCrisprExport() {
  const guideBp = bpFromLength(crGuideLen.value, crGuideUnit.value);
  const insertBp = bpFromLength(crInsertLen.value, crInsertUnit.value);
  const guideConc = parseFloat(crGuideConc.value);
  const insertConc = parseFloat(crInsertConc.value);
  const guideMaxVol = parseFloat(crGuideVol.value);
  const insertMaxVol = parseFloat(crInsertVol.value);
  const guideParts = parseFloat(crGuideParts.value);
  const insertParts = parseFloat(crInsertParts.value);

  const lines = [
    "CRISPR/Cas9 Electroporation Mix Calculator",
    `Generated: ${new Date().toLocaleString()}`,
    "",
    `Guide:  ${crGuideLen.value} ${crGuideUnit.value}, ${crGuideConc.value} ng/µL, ${crGuideVol.value} µL available`,
    `Insert: ${crInsertLen.value} ${crInsertUnit.value}, ${crInsertConc.value} ng/µL, ${crInsertVol.value} µL available`,
    `Ratio (guide:insert): ${crGuideParts.value} : ${crInsertParts.value}`,
    "",
  ];

  const valid = guideBp && insertBp && guideConc > 0 && insertConc > 0 &&
    guideMaxVol > 0 && insertMaxVol > 0 && guideParts > 0 && insertParts > 0;

  if (!valid) {
    lines.push("Result: (enter valid values)");
    return lines.join("\n");
  }

  const guideMaxPmol = pmolFromNg(guideConc * guideMaxVol, guideBp);
  const insertMaxPmol = pmolFromNg(insertConc * insertMaxVol, insertBp);
  const guideUnitPmol = guideMaxPmol / guideParts;
  const insertUnitPmol = insertMaxPmol / insertParts;
  const unitPmol = Math.min(guideUnitPmol, insertUnitPmol);
  const limiting = guideUnitPmol <= insertUnitPmol ? "Guide" : "Insert";
  const guidePmolUsed = unitPmol * guideParts;
  const insertPmolUsed = unitPmol * insertParts;
  const guideNgUsed = ngFromPmol(guidePmolUsed, guideBp);
  const insertNgUsed = ngFromPmol(insertPmolUsed, insertBp);
  const guideVolUsed = guideNgUsed / guideConc;
  const insertVolUsed = insertNgUsed / insertConc;

  lines.push(`Guide to add:  ${fmt(guideVolUsed)} µL  (${fmt(guideNgUsed)} ng · ${fmt(guidePmolUsed, 4)} pmol)`);
  lines.push(`Insert to add: ${fmt(insertVolUsed)} µL  (${fmt(insertNgUsed)} ng · ${fmt(insertPmolUsed, 4)} pmol)`);
  lines.push("", `Limiting factor: ${limiting}`);
  lines.push(`${limiting} uses its full available volume (${fmt(limiting === "Guide" ? guideMaxVol : insertMaxVol)} µL).`);
  lines.push(`Combined DNA volume: ${fmt(guideVolUsed + insertVolUsed)} µL`);

  return lines.join("\n");
}

function buildWesternExport() {
  const lines = [
    "Western Blot Sample Buffer Calculator",
    `Generated: ${new Date().toLocaleString()}`,
    "",
    `Sample buffer factor: ${wbBufferFactor.value}x`,
    `DTT stock factor: ${wbDttFactor.value}x`,
    `Known: ${wbKnownTitle.textContent} = ${wbKnownVol.value} µL`,
    "",
  ];

  const bufferFactor = parseFloat(wbBufferFactor.value);
  const dttFactor = parseFloat(wbDttFactor.value);
  const knownVol = parseFloat(wbKnownVol.value);
  const validFactors = bufferFactor > 1 && dttFactor > 1;
  const sampleFraction = validFactors ? 1 - 1 / bufferFactor - 1 / dttFactor : NaN;

  if (!validFactors || sampleFraction <= 0 || !knownVol || knownVol <= 0) {
    lines.push("Result: (enter valid values)");
    return lines.join("\n");
  }

  const total = wbKnow === "total" ? knownVol : knownVol / sampleFraction;
  const buffer = total / bufferFactor;
  const dtt = total / dttFactor;
  const sample = total - buffer - dtt;

  lines.push(`Sample:           ${fmt(sample)} µL`);
  lines.push(`DTT:              ${fmt(dtt)} µL`);
  lines.push(`4x Sample buffer: ${fmt(buffer)} µL`);
  lines.push(`Total:            ${fmt(total)} µL`);

  return lines.join("\n");
}

function buildPrimerExport() {
  const data = POLY_DATA[prEnzyme.value];
  const lines = [
    "Primer Tm / Ta Calculator",
    `Generated: ${new Date().toLocaleString()}`,
    "",
    `Enzyme: ${data.label}`,
    `Primer concentration: ${prConc.value} nM`,
  ];
  if (prEnzyme.value === "custom") lines.push(`Monovalent salt: ${prSalt.value} mM`);

  const seq1raw = prSeq1.value.trim();
  const seq2raw = prSeq2.value.trim();
  lines.push("", `Primer 1: ${seq1raw || "(none)"}`);
  lines.push(`  Length ${prLen1.textContent}, GC ${prGc1.textContent}, Tm ${prTm1.textContent}`);
  lines.push(`Primer 2: ${seq2raw || "(none)"}`);
  lines.push(`  Length ${prLen2.textContent}, GC ${prGc2.textContent}, Tm ${prTm2.textContent}`);
  lines.push("", `Annealing temp (Ta): ${prTa.textContent}`);
  if (prStatus.textContent) lines.push(`Note: ${prStatus.textContent}`);

  lines.push("", "Dimer Check");
  const appendDimer = (title, container) => {
    lines.push("", title);
    const text = (container.innerText || container.textContent).trim();
    lines.push(...text.split("\n").map((l) => "  " + l));
  };
  appendDimer("Self-Dimer — Primer 1", prDimerSelf1);
  appendDimer("Self-Dimer — Primer 2", prDimerSelf2);
  appendDimer("Hetero-Dimer — Primer 1 × Primer 2", prDimerHetero);

  return lines.join("\n");
}

// ================= RESTRICTION DIGEST =================
// Top-level Test/Mass-based mode switch, scoped to #tool-digest so it can't collide with the
// HiFi tool's own .mode-btn/.mode-panel switch elsewhere on the page.
document.getElementById("digest-mode-switch").querySelectorAll(".mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    document.getElementById("digest-mode-switch").querySelectorAll(".mode-btn").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll("#tool-digest .mode-panel").forEach((p) => p.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById(`digestmode-${btn.dataset.mode}`).classList.add("active");
  });
});

// ---- Test Digest ----
const rdTestEnzcountSwitch = document.getElementById("rd-test-enzcount-switch");
const rdTestTotal = document.getElementById("rd-test-total");
const rdTestDnaVol = document.getElementById("rd-test-dna-vol");
const rdTestEnz1Vol = document.getElementById("rd-test-enz1-vol");
const rdTestEnz2Vol = document.getElementById("rd-test-enz2-vol");
const rdTestEnz2Wrap = document.getElementById("rd-test-enz2-wrap");
const rdTestOutEnz2Row = document.getElementById("rd-test-out-enz2-row");
const rdTestOutDna = document.getElementById("rd-test-out-dna");
const rdTestOutEnz1 = document.getElementById("rd-test-out-enz1");
const rdTestOutEnz2 = document.getElementById("rd-test-out-enz2");
const rdTestOutBuffer = document.getElementById("rd-test-out-buffer");
const rdTestOutWater = document.getElementById("rd-test-out-water");
const rdTestOutTotal = document.getElementById("rd-test-out-total");
const rdTestStatus = document.getElementById("rd-test-status");

let rdTestEnzCount = 1;

rdTestEnzcountSwitch.querySelectorAll(".mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    rdTestEnzcountSwitch.querySelectorAll(".mode-btn").forEach((b) => b.classList.remove("active"));
    btn.classList.add("active");
    rdTestEnzCount = parseInt(btn.dataset.count, 10);
    const show2 = rdTestEnzCount === 2;
    rdTestEnz2Wrap.style.display = show2 ? "" : "none";
    rdTestOutEnz2Row.style.display = show2 ? "" : "none";
    // lab convention: 1 enzyme -> 30 µL rxn / 0.3 µL enzyme; 2 enzymes -> 50 µL rxn / 0.5 µL each
    rdTestTotal.value = show2 ? 50 : 30;
    rdTestEnz1Vol.value = show2 ? 0.5 : 0.3;
    if (show2) rdTestEnz2Vol.value = 0.5;
    calcRdTest();
  });
});

// ---- Test Digest master mix (multiple constructs, same enzymes) ----
const rdTestMmN = document.getElementById("rd-test-mm-n");
const rdTestMmExtra = document.getElementById("rd-test-mm-extra");
const rdTestMmEnz1 = document.getElementById("rd-test-mm-enz1");
const rdTestMmEnz2Row = document.getElementById("rd-test-mm-enz2-row");
const rdTestMmEnz2 = document.getElementById("rd-test-mm-enz2");
const rdTestMmBuffer = document.getElementById("rd-test-mm-buffer");
const rdTestMmWater = document.getElementById("rd-test-mm-water");
const rdTestMmTotal = document.getElementById("rd-test-mm-total");
const rdTestMmAliquot = document.getElementById("rd-test-mm-aliquot");
const rdTestMmStatus = document.getElementById("rd-test-mm-status");

function calcRdTestMasterMix() {
  const total = parseFloat(rdTestTotal.value);
  const dna = parseFloat(rdTestDnaVol.value);
  const enz1 = parseFloat(rdTestEnz1Vol.value);
  const enz2 = rdTestEnzCount === 2 ? parseFloat(rdTestEnz2Vol.value) : 0;
  const n = parseInt(rdTestMmN.value, 10);
  const extra = parseInt(rdTestMmExtra.value, 10) || 0;

  rdTestMmEnz2Row.style.display = rdTestEnzCount === 2 ? "" : "none";

  const clear = () => {
    [rdTestMmEnz1, rdTestMmEnz2, rdTestMmBuffer, rdTestMmWater, rdTestMmTotal, rdTestMmAliquot].forEach((el) => (el.textContent = "—"));
  };

  const valid = total > 0 && dna >= 0 && enz1 >= 0 && n >= 1 && extra >= 0 && (rdTestEnzCount === 1 || enz2 >= 0);
  if (!valid) {
    clear();
    rdTestMmStatus.textContent = "";
    return;
  }

  const buffer = total / 10;
  const enzSum = enz1 + enz2;
  const water = total - dna - enzSum - buffer;

  if (water < 0) {
    clear();
    rdTestMmStatus.textContent = "DNA + enzyme(s) + buffer exceed the reaction volume — fix that above first.";
    return;
  }

  const aliquot = total - dna; // per-reaction master mix volume: everything except each construct's own DNA
  const scale = n + extra;

  rdTestMmEnz1.textContent = `${fmt(enz1 * scale)} µL`;
  if (rdTestEnzCount === 2) rdTestMmEnz2.textContent = `${fmt(enz2 * scale)} µL`;
  rdTestMmBuffer.textContent = `${fmt(buffer * scale)} µL`;
  rdTestMmWater.textContent = `${fmt(water * scale)} µL`;
  rdTestMmTotal.textContent = `${fmt(aliquot * scale)} µL`;
  rdTestMmAliquot.textContent = `${fmt(aliquot)} µL`;

  rdTestMmStatus.textContent =
    `Mixed for ${n} construct${n === 1 ? "" : "s"} + ${extra} extra (${scale} reactions total). ` +
    `Aliquot ${fmt(aliquot)} µL into each of ${n} tubes, then add ${fmt(dna)} µL DNA per construct ` +
    `for a ${fmt(total)} µL final reaction.`;
}

[rdTestMmN, rdTestMmExtra].forEach((el) => el.addEventListener("input", calcRdTestMasterMix));

function calcRdTest() {
  const total = parseFloat(rdTestTotal.value);
  const dna = parseFloat(rdTestDnaVol.value);
  const enz1 = parseFloat(rdTestEnz1Vol.value);
  const enz2 = rdTestEnzCount === 2 ? parseFloat(rdTestEnz2Vol.value) : 0;

  const valid = total > 0 && dna >= 0 && enz1 >= 0 && (rdTestEnzCount === 1 || enz2 >= 0);
  if (!valid) {
    [rdTestOutDna, rdTestOutEnz1, rdTestOutEnz2, rdTestOutBuffer, rdTestOutWater, rdTestOutTotal].forEach((el) => (el.textContent = "—"));
    rdTestStatus.textContent = "";
    calcRdTestMasterMix();
    return;
  }

  const buffer = total / 10;
  const enzSum = enz1 + enz2;
  const water = total - dna - enzSum - buffer;

  rdTestOutDna.textContent = `${fmt(dna)} µL`;
  rdTestOutEnz1.textContent = `${fmt(enz1)} µL`;
  if (rdTestEnzCount === 2) rdTestOutEnz2.textContent = `${fmt(enz2)} µL`;
  rdTestOutBuffer.textContent = `${fmt(buffer)} µL`;
  rdTestOutWater.textContent = `${fmt(water)} µL`;
  rdTestOutTotal.textContent = `${fmt(total)} µL`;

  let status = "";
  if (water < 0) {
    status = "DNA + enzyme(s) + buffer exceed the total reaction volume — increase total volume.";
  } else if (enzSum > 0.1 * total) {
    status = `Enzyme volume is ${fmt((enzSum / total) * 100, 1)}% of the reaction — NEB recommends keeping enzyme ≤10% of total volume to avoid glycerol effects.`;
  }
  rdTestStatus.textContent = status;

  calcRdTestMasterMix();
}

[rdTestTotal, rdTestDnaVol, rdTestEnz1Vol, rdTestEnz2Vol].forEach((el) => el.addEventListener("input", calcRdTest));
calcRdTest();

// ---- Mass-Based Digest ----
const rdMassEnzcountSwitch = document.getElementById("rd-mass-enzcount-switch");
const rdMassUg = document.getElementById("rd-mass-ug");
const rdMassDnaConc = document.getElementById("rd-mass-dna-conc");
const rdMassTotal = document.getElementById("rd-mass-total");
const rdMassUnitsPerUg = document.getElementById("rd-mass-units-per-ug");
const rdMassEnz1Conc = document.getElementById("rd-mass-enz1-conc");
const rdMassEnz2Conc = document.getElementById("rd-mass-enz2-conc");
const rdMassEnz2Wrap = document.getElementById("rd-mass-enz2-wrap");
const rdMassOutEnz2Row = document.getElementById("rd-mass-out-enz2-row");
const rdMassOutDna = document.getElementById("rd-mass-out-dna");
const rdMassOutEnz1 = document.getElementById("rd-mass-out-enz1");
const rdMassOutEnz2 = document.getElementById("rd-mass-out-enz2");
const rdMassOutBuffer = document.getElementById("rd-mass-out-buffer");
const rdMassOutWater = document.getElementById("rd-mass-out-water");
const rdMassOutTotal = document.getElementById("rd-mass-out-total");
const rdMassStatus = document.getElementById("rd-mass-status");

let rdMassEnzCount = 1;

rdMassEnzcountSwitch.querySelectorAll(".mode-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    rdMassEnzcountSwitch.querySelectorAll(".mode-btn").forEach((b) => b.classList.remove("active"));
    btn.classList.add("active");
    rdMassEnzCount = parseInt(btn.dataset.count, 10);
    const show2 = rdMassEnzCount === 2;
    rdMassEnz2Wrap.style.display = show2 ? "" : "none";
    rdMassOutEnz2Row.style.display = show2 ? "" : "none";
    calcRdMass();
  });
});

document.querySelectorAll("#rd-mass-presets .preset-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    rdMassUg.value = btn.dataset.preset;
    calcRdMass();
  });
});

function calcRdMass() {
  const targetUg = parseFloat(rdMassUg.value);
  const dnaConc = parseFloat(rdMassDnaConc.value);
  const total = parseFloat(rdMassTotal.value);
  const unitsPerUg = parseFloat(rdMassUnitsPerUg.value);
  const enz1Conc = parseFloat(rdMassEnz1Conc.value);
  const enz2Conc = rdMassEnzCount === 2 ? parseFloat(rdMassEnz2Conc.value) : null;

  const valid = targetUg > 0 && dnaConc > 0 && total > 0 && unitsPerUg > 0 && enz1Conc > 0 && (rdMassEnzCount === 1 || enz2Conc > 0);
  if (!valid) {
    [rdMassOutDna, rdMassOutEnz1, rdMassOutEnz2, rdMassOutBuffer, rdMassOutWater, rdMassOutTotal].forEach((el) => (el.textContent = "—"));
    rdMassStatus.textContent = "";
    return;
  }

  const dnaVol = (targetUg * 1000) / dnaConc;
  const unitsNeeded = targetUg * unitsPerUg;
  const enz1Vol = unitsNeeded / enz1Conc;
  const enz2Vol = rdMassEnzCount === 2 ? unitsNeeded / enz2Conc : 0;
  const buffer = total / 10;
  const enzSum = enz1Vol + enz2Vol;
  const water = total - dnaVol - enzSum - buffer;

  const lowVolNote = (vol, what) => (vol < 1 ? `<span class="mm-note">⚠ below 1 µL — consider diluting ${what}</span>` : "");

  rdMassOutDna.innerHTML = `${fmt(dnaVol)} µL${lowVolNote(dnaVol, "DNA")}`;
  rdMassOutEnz1.innerHTML = `${fmt(enz1Vol)} µL${lowVolNote(enz1Vol, "enzyme")}`;
  if (rdMassEnzCount === 2) rdMassOutEnz2.innerHTML = `${fmt(enz2Vol)} µL${lowVolNote(enz2Vol, "enzyme")}`;
  rdMassOutBuffer.textContent = `${fmt(buffer)} µL`;
  rdMassOutWater.textContent = `${fmt(water)} µL`;
  rdMassOutTotal.textContent = `${fmt(total)} µL`;

  let status = "";
  if (water < 0) {
    status = "DNA + enzyme(s) + buffer exceed the total reaction volume — increase total volume or use more concentrated stocks.";
  } else if (enzSum > 0.1 * total) {
    status = `Enzyme volume is ${fmt((enzSum / total) * 100, 1)}% of the reaction — NEB recommends keeping enzyme ≤10% of total volume to avoid glycerol effects.`;
  }
  rdMassStatus.textContent = status;
}

[rdMassUg, rdMassDnaConc, rdMassTotal, rdMassUnitsPerUg, rdMassEnz1Conc, rdMassEnz2Conc].forEach((el) =>
  el.addEventListener("input", calcRdMass)
);
calcRdMass();

function buildDigestExport() {
  const isTest = document.getElementById("digestmode-test").classList.contains("active");
  const lines = ["Restriction Digest Calculator", `Generated: ${new Date().toLocaleString()}`, ""];

  if (isTest) {
    lines.push(`Mode: Test Digest (Diagnostic) — ${rdTestEnzCount} enzyme(s)`, "");
    const total = parseFloat(rdTestTotal.value);
    const dna = parseFloat(rdTestDnaVol.value);
    const enz1 = parseFloat(rdTestEnz1Vol.value);
    const enz2 = rdTestEnzCount === 2 ? parseFloat(rdTestEnz2Vol.value) : 0;

    if (!(total > 0 && dna >= 0 && enz1 >= 0)) {
      lines.push("Result: (enter valid values)");
      return lines.join("\n");
    }

    const buffer = total / 10;
    const water = total - dna - enz1 - enz2 - buffer;

    lines.push(`DNA: ${fmt(dna)} µL`);
    lines.push(`Enzyme 1: ${fmt(enz1)} µL`);
    if (rdTestEnzCount === 2) lines.push(`Enzyme 2: ${fmt(enz2)} µL`);
    lines.push(`10x NEB Buffer: ${fmt(buffer)} µL`);
    lines.push(`ddH2O: ${fmt(water)} µL`);
    lines.push(`Total: ${fmt(total)} µL`);
    if (rdTestStatus.textContent) lines.push("", `Note: ${rdTestStatus.textContent}`);

    const n = parseInt(rdTestMmN.value, 10);
    const extra = parseInt(rdTestMmExtra.value, 10) || 0;
    if (n >= 1 && water >= 0) {
      const aliquot = total - dna;
      const scale = n + extra;
      lines.push("", `Master Mix — ${n} construct(s) + ${extra} extra (${scale} reactions total)`);
      lines.push(`  Enzyme 1: ${fmt(enz1 * scale)} µL`);
      if (rdTestEnzCount === 2) lines.push(`  Enzyme 2: ${fmt(enz2 * scale)} µL`);
      lines.push(`  10x NEB Buffer: ${fmt(buffer * scale)} µL`);
      lines.push(`  ddH2O: ${fmt(water * scale)} µL`);
      lines.push(`  Master mix total: ${fmt(aliquot * scale)} µL`);
      lines.push(`  Aliquot per tube: ${fmt(aliquot)} µL, then add ${fmt(dna)} µL DNA per construct`);
    }
  } else {
    lines.push(`Mode: Mass-Based Digest — ${rdMassEnzCount} enzyme(s)`, "");
    const targetUg = parseFloat(rdMassUg.value);
    const dnaConc = parseFloat(rdMassDnaConc.value);
    const total = parseFloat(rdMassTotal.value);
    const unitsPerUg = parseFloat(rdMassUnitsPerUg.value);
    const enz1Conc = parseFloat(rdMassEnz1Conc.value);
    const enz2Conc = rdMassEnzCount === 2 ? parseFloat(rdMassEnz2Conc.value) : null;

    if (!(targetUg > 0 && dnaConc > 0 && total > 0 && unitsPerUg > 0 && enz1Conc > 0 && (rdMassEnzCount === 1 || enz2Conc > 0))) {
      lines.push("Result: (enter valid values)");
      return lines.join("\n");
    }

    const dnaVol = (targetUg * 1000) / dnaConc;
    const unitsNeeded = targetUg * unitsPerUg;
    const enz1Vol = unitsNeeded / enz1Conc;
    const enz2Vol = rdMassEnzCount === 2 ? unitsNeeded / enz2Conc : 0;
    const buffer = total / 10;
    const water = total - dnaVol - enz1Vol - enz2Vol - buffer;

    lines.push(`Target DNA: ${targetUg} µg at ${dnaConc} ng/µL → ${fmt(dnaVol)} µL`);
    lines.push(`Units needed per enzyme: ${fmt(unitsNeeded)} U`);
    lines.push(`Enzyme 1 (${enz1Conc} U/µL): ${fmt(enz1Vol)} µL`);
    if (rdMassEnzCount === 2) lines.push(`Enzyme 2 (${enz2Conc} U/µL): ${fmt(enz2Vol)} µL`);
    lines.push(`10x NEB Buffer: ${fmt(buffer)} µL`);
    lines.push(`ddH2O: ${fmt(water)} µL`);
    lines.push(`Total: ${fmt(total)} µL`);
    if (rdMassStatus.textContent) lines.push("", `Note: ${rdMassStatus.textContent}`);
  }

  return lines.join("\n");
}

const EXPORT_BUILDERS = {
  hifi: buildHifiExport,
  crispr: buildCrisprExport,
  western: buildWesternExport,
  primer: buildPrimerExport,
  digest: buildDigestExport,
  precip: buildPrecipExport,
};

document.querySelectorAll(".export-btn").forEach((btn) => {
  btn.addEventListener("click", () => {
    const tool = btn.dataset.tool;
    const text = EXPORT_BUILDERS[tool]();
    downloadText(exportFilename(tool), text);
  });
});

// ================= SODIUM ACETATE / ETHANOL PRECIPITATION =================
const pcSampleVol = document.getElementById("pc-sample-vol");
const pcOutSample = document.getElementById("pc-out-sample");
const pcOutNaoac = document.getElementById("pc-out-naoac");
const pcOutPretotal = document.getElementById("pc-out-pretotal");
const pcOutEtoh = document.getElementById("pc-out-etoh");
const pcOutTotal = document.getElementById("pc-out-total");

function calcPrecip() {
  const sample = parseFloat(pcSampleVol.value);

  if (!(sample > 0)) {
    [pcOutSample, pcOutNaoac, pcOutPretotal, pcOutEtoh, pcOutTotal].forEach((el) => (el.textContent = "—"));
    return;
  }

  const naoac = sample / 10;
  const pretotal = sample + naoac;
  const etoh = sample * 2.5;
  const total = pretotal + etoh;

  pcOutSample.textContent = `${fmt(sample)} µL`;
  pcOutNaoac.textContent = `${fmt(naoac)} µL`;
  pcOutPretotal.textContent = `${fmt(pretotal)} µL`;
  pcOutEtoh.textContent = `${fmt(etoh)} µL`;
  pcOutTotal.textContent = `${fmt(total)} µL`;
}

pcSampleVol.addEventListener("input", calcPrecip);
calcPrecip();

function buildPrecipExport() {
  const lines = [
    "Sodium Acetate / Ethanol Precipitation Calculator",
    `Generated: ${new Date().toLocaleString()}`,
    "",
  ];

  const sample = parseFloat(pcSampleVol.value);
  if (!(sample > 0)) {
    lines.push("Result: (enter a valid sample volume)");
    return lines.join("\n");
  }

  const naoac = sample / 10;
  const pretotal = sample + naoac;
  const etoh = sample * 2.5;
  const total = pretotal + etoh;

  lines.push(`Sample: ${fmt(sample)} µL`);
  lines.push(`3M sodium acetate, pH 5.2 (1/10 vol): ${fmt(naoac)} µL`);
  lines.push(`Volume before ethanol: ${fmt(pretotal)} µL`);
  lines.push(`Ice-cold 100% ethanol (2.5x sample vol): ${fmt(etoh)} µL`);
  lines.push(`Total volume in tube: ${fmt(total)} µL`);

  lines.push("", "Protocol:");
  lines.push("  1. Add 3M sodium acetate, pH 5.2 (1/10th sample volume) to the sample.");
  lines.push("  2. Add ice-cold 100% ethanol (2.5x sample volume) and mix well by inverting.");
  lines.push("  3. Put samples at -20C overnight.");
  lines.push("  4. Spin at 21,000 x g, 4C, for 30 min.");
  lines.push("  5. Decant supernatant by pouring - it's okay if some ethanol is left at the bottom.");
  lines.push("  6. Wash pellet in 500 uL ice-cold 75% ethanol.");
  lines.push("  7. Spin at 21,000 x g, 4C, for 10 min.");
  lines.push("  8. Decant supernatant and dry pellet for 10-15 min.");
  lines.push("  9. Resuspend pellet in water once slightly damp. If overdried, heat at 55C to solubilize - make sure it is all solubilized.");

  return lines.join("\n");
}

// ================= GUIDE DESIGN (links out) =================
const KOGUIDE_URL = "http://127.0.0.1:5001/";
const EUPAGDT_URL = "http://grna.ctegd.uga.edu/";

document.getElementById("koguide-open").addEventListener("click", () => {
  window.open(KOGUIDE_URL, "_blank", "noopener");
});
document.getElementById("eupagdt-open").addEventListener("click", () => {
  window.open(EUPAGDT_URL, "_blank", "noopener");
});

async function checkKoguideStatus() {
  const statusEl = document.getElementById("koguide-status");
  try {
    await fetch(KOGUIDE_URL, { mode: "no-cors", cache: "no-store", signal: AbortSignal.timeout(1500) });
    statusEl.textContent = "● Running — the button below will open it";
    statusEl.classList.add("status-ok");
  } catch (e) {
    statusEl.textContent = '○ Not running — start "KO Guide Design" from the Desktop first';
    statusEl.classList.remove("status-ok");
  }
}

checkKoguideStatus();
