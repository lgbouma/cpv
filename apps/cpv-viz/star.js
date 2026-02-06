const PDF_MANIFEST_URL = "./data/cpv_pdf_manifest.json";
const DATA_MANIFEST_URL = "./data/cpv_lit_subset_manifest.json";
const FALLBACK_DATA_URL = "./data/cpv_lit_subset.json";
const FALLBACK_BASE_URL = "./pdfs";

const dom = {
  starTitle: document.getElementById("star-title"),
  starMeta: document.getElementById("star-meta"),
  pdfSubtitle: document.getElementById("pdf-subtitle"),
  pdfCount: document.getElementById("pdf-count"),
  pdfList: document.getElementById("pdf-list"),
  pdfEmpty: document.getElementById("pdf-empty"),
  detailGrid: document.getElementById("detail-grid"),
  detailEmpty: document.getElementById("detail-empty"),
};

const numberFormatter = new Intl.NumberFormat("en-US", {
  maximumFractionDigits: 3,
});

const plainNumberKeys = new Set(["ticid"]);

function getTicidFromQuery() {
  const params = new URLSearchParams(window.location.search);
  const ticid = params.get("ticid");
  if (!ticid) {
    return null;
  }
  return ticid.trim();
}

function formatValue(value, key) {
  if (value === null || value === undefined || value === "") {
    return "—";
  }
  if (key && plainNumberKeys.has(key)) {
    return String(value);
  }
  if (typeof value === "number") {
    return numberFormatter.format(value);
  }
  return String(value);
}

function extractSectorId(filename) {
  const parts = filename.split("_");
  if (parts.length < 2) {
    return null;
  }
  const match = parts[1].match(/\d+/g);
  if (!match) {
    return null;
  }
  return Number(match.join(""));
}

function sortPdfs(pdfs) {
  return [...pdfs].sort((a, b) => {
    const aSector = extractSectorId(a);
    const bSector = extractSectorId(b);
    if (aSector !== null && bSector !== null && aSector !== bSector) {
      return aSector - bSector;
    }
    if (aSector !== null && bSector === null) {
      return -1;
    }
    if (aSector === null && bSector !== null) {
      return 1;
    }
  return a.localeCompare(b);
  });
}

function updateHeader(ticid, pdfCount) {
  dom.starTitle.textContent = ticid ? `TIC ${ticid} PDFs` : "CPV Vetting PDFs";
  dom.starMeta.textContent = ticid
    ? "Scroll to review the available CPV vetter PDFs for this star."
    : "Provide a TIC ID to view CPV vetter PDFs.";
  dom.pdfCount.textContent = `${pdfCount} PDFs`;
}

function updateSubtitle(manifest) {
  if (!manifest) {
    dom.pdfSubtitle.textContent = "Sorted by sector.";
    return;
  }
  const generated = manifest.generated_date ? `Manifest updated ${manifest.generated_date}.` : null;
  dom.pdfSubtitle.textContent = generated ? `Sorted by sector. ${generated}` : "Sorted by sector.";
}

function normalizeBaseUrl(baseUrl) {
  if (!baseUrl) {
    return FALLBACK_BASE_URL;
  }
  return baseUrl.endsWith("/") ? baseUrl.slice(0, -1) : baseUrl;
}

function renderDetails(row, pdfCount) {
  dom.detailGrid.innerHTML = "";

  if (!row) {
    dom.detailEmpty.hidden = false;
    return;
  }

  dom.detailEmpty.hidden = true;
  const fields = [
    { key: "original_id", label: "Original ID", value: row.original_id },
    { key: "ticid", label: "TIC ID", value: row.ticid },
    { key: "pdf_count", label: "PDFs", value: pdfCount },
    { key: "N_sectors", label: "N_sectors", value: row.N_sectors },
    { key: "distance_pc", label: "Distance (pc)", value: row.distance_pc },
    { key: "TESSMAG", label: "TESSMag", value: row.TESSMAG },
    { key: "cluster", label: "Cluster", value: row.cluster },
    { key: "bibcode", label: "Bibcode", value: row.bibcode },
    { key: "period_hr", label: "Period (hr)", value: row.period_hr },
    { key: "quality", label: "Quality", value: row.quality },
    { key: "telescope", label: "Telescope", value: row.telescope },
  ];

  const fragment = document.createDocumentFragment();
  fields.forEach((field) => {
    const item = document.createElement("div");
    item.className = "detail-item";

    const label = document.createElement("div");
    label.className = "detail-label";
    label.textContent = field.label;

    const value = document.createElement("div");
    value.className = "detail-value";
    value.textContent = formatValue(field.value, field.key);

    item.appendChild(label);
    item.appendChild(value);
    fragment.appendChild(item);
  });

  dom.detailGrid.appendChild(fragment);
}

function buildPdfCard(filename, baseUrl) {
  const card = document.createElement("article");
  card.className = "pdf-card";

  const header = document.createElement("div");
  header.className = "pdf-card-header";

  const sectorId = extractSectorId(filename);
  const title = document.createElement("h3");
  title.textContent = sectorId !== null ? `Sector ${sectorId}` : "CPV vetter PDF";

  const meta = document.createElement("div");
  meta.className = "pdf-card-meta";
  meta.textContent = filename;

  const link = document.createElement("a");
  link.className = "pdf-link";
  link.href = `${baseUrl}/${encodeURIComponent(filename)}`;
  link.target = "_blank";
  link.rel = "noopener";
  link.textContent = "Open PDF";

  header.appendChild(title);
  header.appendChild(link);

  const iframe = document.createElement("iframe");
  iframe.className = "pdf-frame";
  iframe.src = `${baseUrl}/${encodeURIComponent(filename)}`;
  iframe.title = filename;
  iframe.loading = "lazy";

  card.appendChild(header);
  card.appendChild(meta);
  card.appendChild(iframe);

  return card;
}

async function loadManifest() {
  const response = await fetch(PDF_MANIFEST_URL);
  if (!response.ok) {
    throw new Error(`Failed to load PDF manifest: ${response.status}`);
  }
  return response.json();
}

async function loadDataManifest() {
  let dataUrl = FALLBACK_DATA_URL;

  try {
    const response = await fetch(DATA_MANIFEST_URL);
    if (!response.ok) {
      throw new Error(`Failed to load data manifest: ${response.status}`);
    }
    const manifest = await response.json();
    if (manifest && manifest.data_file) {
      dataUrl = `./data/${manifest.data_file}`;
    }
  } catch (error) {
    dataUrl = FALLBACK_DATA_URL;
  }

  return dataUrl;
}

async function init() {
  const ticid = getTicidFromQuery();
  if (!ticid) {
    updateHeader(null, 0);
    dom.pdfEmpty.hidden = false;
    dom.pdfEmpty.textContent = "Missing TIC ID in the URL. Use ?ticid=123456.";
    return;
  }

  let manifest = null;
  let dataRows = [];
  try {
    manifest = await loadManifest();
  } catch (error) {
    updateHeader(ticid, 0);
    dom.starMeta.textContent = "Unable to load the PDF manifest.";
    dom.pdfEmpty.hidden = false;
    dom.pdfEmpty.textContent =
      "Unable to load data/cpv_pdf_manifest.json. Run build_pdf_manifest.py to generate it.";
    return;
  }

  updateSubtitle(manifest);

  const baseUrl = normalizeBaseUrl(manifest.pdf_base_url || FALLBACK_BASE_URL);
  const pdfs = sortPdfs((manifest.pdfs && manifest.pdfs[ticid]) || []);
  updateHeader(ticid, pdfs.length);

  const dataUrl = await loadDataManifest();
  try {
    const response = await fetch(dataUrl);
    if (!response.ok) {
      throw new Error(`Failed to load data: ${response.status}`);
    }
    dataRows = await response.json();
  } catch (error) {
    dataRows = [];
  }

  const detailRow = dataRows.find((row) => String(row.ticid) === String(ticid)) || null;
  renderDetails(detailRow, pdfs.length);

  dom.pdfList.innerHTML = "";
  if (pdfs.length === 0) {
    dom.pdfEmpty.hidden = false;
    dom.pdfEmpty.textContent = "No PDFs found for this TIC ID.";
    return;
  }

  dom.pdfEmpty.hidden = true;
  const fragment = document.createDocumentFragment();
  pdfs.forEach((filename) => {
    fragment.appendChild(buildPdfCard(filename, baseUrl));
  });
  dom.pdfList.appendChild(fragment);
}

document.addEventListener("DOMContentLoaded", init);
