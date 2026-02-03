const MANIFEST_URL = "./data/cpv_lit_subset_manifest.json";
const FALLBACK_DATA_URL = "./data/cpv_lit_subset.json";

const columnDefs = [
  { key: "original_id", label: "Original ID", type: "string" },
  { key: "ticid", label: "TIC ID", type: "number" },
  { key: "sectors", label: "Sectors", type: "string" },
  { key: "N_sectors", label: "N_sectors", type: "number" },
  { key: "distance_pc", label: "Distance (pc)", type: "number" },
  { key: "TESSMAG", label: "TESSMAG", type: "number" },
  { key: "cluster", label: "Cluster", type: "string" },
  { key: "bibcode", label: "Bibcode", type: "string" },
  { key: "period_hr", label: "Period (hr)", type: "number" },
  { key: "quality", label: "Quality", type: "number" },
  { key: "telescope", label: "Telescope", type: "string" },
];

const state = {
  rows: [],
  sort: { key: "original_id", direction: "asc" },
  filters: {
    cluster: "All",
    quality: "All",
    telescope: "All",
  },
  searchQuery: "",
  showSectors: false,
};

const dom = {
  filterCluster: document.getElementById("filter-cluster"),
  filterQuality: document.getElementById("filter-quality"),
  filterTelescope: document.getElementById("filter-telescope"),
  searchInput: document.getElementById("search-input"),
  toggleSectors: document.getElementById("toggle-sectors"),
  resetFilters: document.getElementById("reset-filters"),
  tableHeaders: document.getElementById("table-headers"),
  tableBody: document.getElementById("table-body"),
  rowCount: document.getElementById("row-count"),
  activeSort: document.getElementById("active-sort"),
  dataUpdated: document.getElementById("data-updated"),
  currentSector: document.getElementById("current-sector"),
};

const numberFormatter = new Intl.NumberFormat("en-US", {
  maximumFractionDigits: 3,
});

const plainNumberKeys = new Set(["ticid"]);

function getVisibleColumns() {
  const baseColumns = columnDefs.filter((column) => column.key !== "sectors");
  if (!state.showSectors) {
    return baseColumns;
  }
  const sectorsColumn = columnDefs.find((column) => column.key === "sectors");
  return sectorsColumn ? [...baseColumns, sectorsColumn] : baseColumns;
}

function formatValue(value, type, key) {
  if (value === null || value === undefined || value === "") {
    return "—";
  }
  if (key && plainNumberKeys.has(key)) {
    return String(value);
  }
  if (type === "number" && typeof value === "number") {
    return numberFormatter.format(value);
  }
  return String(value);
}

function setSort(key) {
  if (state.sort.key === key) {
    state.sort.direction = state.sort.direction === "asc" ? "desc" : "asc";
  } else {
    state.sort.key = key;
    state.sort.direction = "asc";
  }
  render();
}

function buildTableHeaders() {
  dom.tableHeaders.innerHTML = "";
  getVisibleColumns().forEach((column) => {
    const th = document.createElement("th");
    const button = document.createElement("button");
    button.type = "button";
    button.className = "header-button";
    button.dataset.key = column.key;
    button.dataset.label = column.label;
    button.textContent = column.label;
    button.addEventListener("click", () => setSort(column.key));
    th.appendChild(button);
    dom.tableHeaders.appendChild(th);
  });
}

function uniqueValues(key) {
  const values = new Set();
  state.rows.forEach((row) => {
    const value = row[key];
    if (value !== null && value !== undefined && value !== "") {
      values.add(String(value));
    }
  });
  return Array.from(values).sort((a, b) => a.localeCompare(b));
}

function buildFilterButtons(container, key, values) {
  container.innerHTML = "";
  const allButton = document.createElement("button");
  allButton.type = "button";
  allButton.className = "chip-button";
  allButton.textContent = "All";
  allButton.addEventListener("click", () => {
    state.filters[key] = "All";
    render();
  });
  container.appendChild(allButton);

  values.forEach((value) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "chip-button";
    button.textContent = value;
    button.addEventListener("click", () => {
      state.filters[key] = value;
      render();
    });
    container.appendChild(button);
  });
}

function buildFilterSelect(selectEl, key, values) {
  selectEl.innerHTML = "";
  const allOption = document.createElement("option");
  allOption.value = "All";
  allOption.textContent = "All";
  selectEl.appendChild(allOption);

  values.forEach((value) => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value;
    selectEl.appendChild(option);
  });

  selectEl.addEventListener("change", (event) => {
    state.filters[key] = event.target.value;
    render();
  });
}

function updateFilterUI() {
  const buttonGroups = [
    { container: dom.filterQuality, key: "quality" },
    { container: dom.filterTelescope, key: "telescope" },
  ];

  buttonGroups.forEach(({ container, key }) => {
    Array.from(container.querySelectorAll("button")).forEach((button) => {
      const isActive = button.textContent === state.filters[key];
      button.classList.toggle("active", isActive);
    });
  });

  if (dom.filterCluster) {
    dom.filterCluster.value = state.filters.cluster;
  }

  dom.searchInput.value = state.searchQuery;
  dom.toggleSectors.checked = state.showSectors;
}

function updateSortUI() {
  const updateButtons = (buttons) => {
    Array.from(buttons).forEach((button) => {
      const isActive = button.dataset.key === state.sort.key;
      const label = button.dataset.label ?? button.textContent;
      const arrow = isActive ? (state.sort.direction === "asc" ? " ↑" : " ↓") : "";
      button.textContent = `${label}${arrow}`;
      button.classList.toggle("active", isActive);
    });
  };

  updateButtons(dom.tableHeaders.querySelectorAll("button"));

  const activeColumn = columnDefs.find((col) => col.key === state.sort.key);
  const activeLabel = activeColumn ? activeColumn.label : state.sort.key;
  dom.activeSort.textContent = `Sort: ${activeLabel} (${state.sort.direction})`;
}

function updateDataUpdated(displayDate, currentSector) {
  if (!dom.dataUpdated) {
    return;
  }
  const dateText = displayDate || "—";
  dom.dataUpdated.textContent = `Sector counts last updated ${dateText}.`;

  if (dom.currentSector) {
    const sectorText = currentSector ? `S${currentSector}` : "—";
    dom.currentSector.textContent = `Current TESS sector ${sectorText}.`;
  }
}

function applyFilters(rows) {
  const query = state.searchQuery.trim().toLowerCase();
  return rows.filter((row) => {
    if (state.filters.cluster !== "All" && String(row.cluster) !== state.filters.cluster) {
      return false;
    }
    if (state.filters.quality !== "All" && String(row.quality) !== state.filters.quality) {
      return false;
    }
    if (state.filters.telescope !== "All" && String(row.telescope) !== state.filters.telescope) {
      return false;
    }

    if (query) {
      const haystack = columnDefs
        .map((column) => formatValue(row[column.key], column.type, column.key))
        .join(" ")
        .toLowerCase();
      return haystack.includes(query);
    }

    return true;
  });
}

function sortRows(rows) {
  const { key, direction } = state.sort;
  const column = columnDefs.find((col) => col.key === key);
  const multiplier = direction === "asc" ? 1 : -1;
  return [...rows].sort((a, b) => {
    const aVal = a[key];
    const bVal = b[key];

    if (aVal === null || aVal === undefined) {
      return bVal === null || bVal === undefined ? 0 : 1;
    }
    if (bVal === null || bVal === undefined) {
      return -1;
    }

    if (column && column.type === "number") {
      return (aVal - bVal) * multiplier;
    }

    return String(aVal).localeCompare(String(bVal)) * multiplier;
  });
}

function renderTable(rows) {
  const visibleColumns = getVisibleColumns();
  dom.tableBody.innerHTML = "";

  if (rows.length === 0) {
    const emptyRow = document.createElement("tr");
    const emptyCell = document.createElement("td");
    emptyCell.colSpan = visibleColumns.length;
    emptyCell.className = "empty-state";
    emptyCell.textContent = "No rows match the current filters.";
    emptyRow.appendChild(emptyCell);
    dom.tableBody.appendChild(emptyRow);
    return;
  }

  const fragment = document.createDocumentFragment();
  rows.forEach((row) => {
    const tr = document.createElement("tr");
    visibleColumns.forEach((column) => {
      const td = document.createElement("td");
      td.textContent = formatValue(row[column.key], column.type, column.key);
      tr.appendChild(td);
    });
    fragment.appendChild(tr);
  });
  dom.tableBody.appendChild(fragment);
}

function ensureSortKeyVisible() {
  const visibleKeys = new Set(getVisibleColumns().map((column) => column.key));
  if (!visibleKeys.has(state.sort.key)) {
    state.sort.key = "original_id";
    state.sort.direction = "asc";
  }
}

function render() {
  ensureSortKeyVisible();
  const filteredRows = applyFilters(state.rows);
  const sortedRows = sortRows(filteredRows);

  renderTable(sortedRows);
  dom.rowCount.textContent = `${sortedRows.length} rows`;
  updateFilterUI();
  updateSortUI();
}

function bindInputs() {
  dom.searchInput.addEventListener("input", (event) => {
    state.searchQuery = event.target.value;
    render();
  });

  dom.toggleSectors.addEventListener("change", (event) => {
    state.showSectors = event.target.checked;
    buildTableHeaders();
    render();
  });

  dom.resetFilters.addEventListener("click", () => {
    state.filters.cluster = "All";
    state.filters.quality = "All";
    state.filters.telescope = "All";
    state.searchQuery = "";
    render();
  });
}

async function loadManifest() {
  let dataUrl = FALLBACK_DATA_URL;
  let displayDate = null;
  let currentSector = null;

  try {
    const response = await fetch(MANIFEST_URL);
    if (!response.ok) {
      throw new Error(`Failed to load manifest: ${response.status}`);
    }
    const manifest = await response.json();
    if (manifest && manifest.data_file) {
      dataUrl = `./data/${manifest.data_file}`;
    }
    displayDate = manifest.display_date || manifest.generated_date || null;
    currentSector = manifest.current_sector ?? null;
  } catch (error) {
    displayDate = null;
    currentSector = null;
  }

  updateDataUpdated(displayDate, currentSector);
  return dataUrl;
}

async function init() {
  const dataUrl = await loadManifest();

  try {
    const response = await fetch(dataUrl);
    if (!response.ok) {
      throw new Error(`Failed to load data: ${response.status}`);
    }
    state.rows = await response.json();
  } catch (error) {
    dom.tableBody.innerHTML = "";
    const errorRow = document.createElement("tr");
    const errorCell = document.createElement("td");
    errorCell.colSpan = getVisibleColumns().length;
    errorCell.className = "empty-state";
    errorCell.textContent = "Unable to load data. Please check the JSON file.";
    errorRow.appendChild(errorCell);
    dom.tableBody.appendChild(errorRow);
    return;
  }

  buildTableHeaders();
  buildFilterSelect(dom.filterCluster, "cluster", uniqueValues("cluster"));
  buildFilterButtons(dom.filterQuality, "quality", uniqueValues("quality"));
  buildFilterButtons(dom.filterTelescope, "telescope", uniqueValues("telescope"));
  bindInputs();
  render();
}

document.addEventListener("DOMContentLoaded", init);
