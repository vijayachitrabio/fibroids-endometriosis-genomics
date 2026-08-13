const state = {
  manifest: null,
  activeTab: "about",
  tableCache: new Map(),
  search: ""
};

const $ = (selector) => document.querySelector(selector);

const analyses = [
  {
    id: "about",
    label: "About",
    title: "About this results browser",
    text: "This companion browser organizes the manuscript results into static analysis tabs. It is designed for GitHub Pages and contains only aggregate outputs, selected figures, and manuscript-ready tables.",
    bullets: [
      "No UK Biobank participant-level records are included.",
      "No individual phenotype dropdowns or patient-level browsing are provided.",
      "Cancer analyses are intentionally excluded from the main browser tables."
    ],
    figures: ["fig1"],
    tables: []
  },
  {
    id: "phenotype",
    label: "Comorbidity",
    title: "Phenotypic comorbidity architecture",
    text: "Aggregate odds-ratio results summarize how uterine fibroids and endometriosis differ across non-cancer comorbidity domains.",
    figures: ["fig2"],
    tables: ["table2"]
  },
  {
    id: "temporal",
    label: "Temporal",
    title: "Temporal ordering of comorbidities",
    text: "HES diagnosis dates were used to classify comorbidities relative to the primary disease index date using aggregate event counts.",
    figures: ["fig3"],
    tables: ["table3"]
  },
  {
    id: "genetics",
    label: "Genetics",
    title: "Genome-wide and local genetic correlation",
    text: "The genetics tab focuses on LDSC/LAVA evidence for shared local genetic architecture between uterine fibroids and endometriosis.",
    figures: ["fig4"],
    tables: ["s2"]
  },
  {
    id: "magma",
    label: "MAGMA",
    title: "Gene-level architecture",
    text: "MAGMA gene-level results identify disease-specific and shared gene signals, including the 10 shared gene-level associations.",
    figures: ["fig5"],
    tables: ["s14c"]
  },
  {
    id: "mr",
    label: "MR",
    title: "Mendelian randomization and sensitivity analyses",
    text: "This tab shows no-cancer PheWAS-MR results and the MR-PRESSO sensitivity table for feasible multi-instrument pairs.",
    figures: ["fig7"],
    tables: ["s6", "s19"]
  },
  {
    id: "proteomics",
    label: "Proteomics",
    title: "Olink and protein-MR evidence",
    text: "Corrected protein-MR results prioritize TFPI and FLT3LG with hypertension while noting sparse cis-pQTL instrument limitations.",
    figures: ["fig7"],
    tables: ["s8"]
  },
  {
    id: "singlecell",
    label: "Single-cell",
    title: "Single-cell context",
    text: "Public scRNA-seq datasets provide cell-type context for prioritized genes and TFPI expression across uterine fibroid and endometriosis datasets.",
    figures: ["figS"],
    tables: ["s18"]
  },
  {
    id: "downloads",
    label: "Downloads",
    title: "Download aggregate tables",
    text: "All bundled browser tables are listed here for direct retrieval. These files are static aggregate outputs only.",
    figures: [],
    tables: "all"
  }
];

function el(tag, options = {}, children = []) {
  const node = document.createElement(tag);
  Object.entries(options).forEach(([key, value]) => {
    if (key === "className") node.className = value;
    else if (key === "text") node.textContent = value;
    else if (key === "html") node.innerHTML = value;
    else node.setAttribute(key, value);
  });
  children.forEach((child) => node.appendChild(child));
  return node;
}

function parseDelimited(text, delimiter) {
  const rows = [];
  let row = [];
  let cell = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    const next = text[i + 1];
    if (char === '"' && quoted && next === '"') {
      cell += '"';
      i += 1;
    } else if (char === '"') {
      quoted = !quoted;
    } else if (char === delimiter && !quoted) {
      row.push(cell);
      cell = "";
    } else if ((char === "\n" || char === "\r") && !quoted) {
      if (char === "\r" && next === "\n") i += 1;
      row.push(cell);
      if (row.some((value) => value.length > 0)) rows.push(row);
      row = [];
      cell = "";
    } else {
      cell += char;
    }
  }
  row.push(cell);
  if (row.some((value) => value.length > 0)) rows.push(row);
  if (!rows.length) return { columns: [], records: [] };
  const columns = rows[0].map((value) => value.trim());
  const records = rows.slice(1).map((values) => {
    const record = {};
    columns.forEach((column, index) => {
      record[column] = values[index] ?? "";
    });
    return record;
  });
  return { columns, records };
}

function tableById(id) {
  return state.manifest.tables.find((table) => table.id === id);
}

function figureById(id) {
  return state.manifest.figures.find((figure) => figure.id === id);
}

function compact(value) {
  const text = String(value ?? "");
  if (text.length <= 90) return text;
  return `${text.slice(0, 87)}...`;
}

async function loadTableData(table) {
  if (state.tableCache.has(table.id)) return state.tableCache.get(table.id);
  const response = await fetch(table.path);
  if (!response.ok) throw new Error(`Could not load ${table.path}`);
  const text = await response.text();
  const parsed = parseDelimited(text, table.format === "tsv" ? "\t" : ",");
  state.tableCache.set(table.id, parsed);
  return parsed;
}

function renderTabs() {
  const tabs = $("#tabs");
  tabs.innerHTML = "";
  analyses.forEach((analysis) => {
    const li = el("li", { className: "nav-item", role: "presentation" });
    
    // Create button matching the user's screenshot (white text, active has translucent background)
    const isActive = analysis.id === state.activeTab;
    const button = el("button", {
      className: `btn fw-medium px-3 py-2 text-white ${isActive ? "bg-white bg-opacity-25" : "bg-transparent"} border-0 rounded-2`,
      type: "button",
      role: "tab",
      "aria-selected": String(isActive),
      text: analysis.label
    });
    
    button.addEventListener("click", () => {
      state.activeTab = analysis.id;
      state.search = "";
      render();
      $("#app").focus();
    });
    
    li.appendChild(button);
    tabs.appendChild(li);
  });
}

function renderMetrics() {
  const list = el("div", { className: "row g-3 mt-4" });
  state.manifest.metrics.forEach((metric) => {
    list.appendChild(el("div", { className: "col-sm-6 col-lg-4" }, [
      el("div", { className: "card h-100 bg-light border-0 shadow-sm" }, [
        el("div", { className: "card-body p-3" }, [
          el("h3", { className: "card-title text-primary fw-bold mb-1 fs-4", text: metric.value }),
          el("p", { className: "card-text small text-muted mb-0", text: metric.label })
        ])
      ])
    ]));
  });
  return list;
}

function renderAbout(analysis) {
  const blocks = [
    ["Included", "Curated aggregate figures, static CSV/TSV tables, no-cancer MR summaries, MR-PRESSO sensitivity results, and public scRNA-seq summaries."],
    ["Excluded", "UK Biobank participant-level data, raw phenotype matrices, patient-level browsers, and exploratory cancer analysis tabs."],
    ["How to use", "Select an analysis tab, inspect the manuscript figure, search the accompanying table, or download the static aggregate file."]
  ];
  
  const aboutLayout = el("div", { className: "row g-4 mb-4" }, [
    el("div", { className: "col-lg-8" }, [
      el("div", { className: "card shadow-sm h-100" }, [
        el("div", { className: "card-body p-4" }, [
          el("span", { className: "badge bg-info text-dark mb-2", text: "About" }),
          el("h2", { className: "card-title fw-bold mb-3", text: analysis.title }),
          el("p", { className: "card-text lead text-muted", text: analysis.text }),
          el("div", { className: "alert alert-warning mt-4", role: "alert", text: state.manifest.complianceNote }),
          renderMetrics()
        ])
      ])
    ]),
    el("div", { className: "col-lg-4" }, [
      el("div", { className: "card shadow-sm bg-light border-0 h-100" }, [
        el("div", { className: "card-body p-4" }, [
          el("h3", { className: "h5 fw-bold text-dark border-bottom pb-2 mb-3", text: "Use boundary" }),
          el("ul", { className: "ps-3 text-muted mb-4" }, analysis.bullets.map((item) => el("li", { className: "mb-2", text: item }))),
          ...blocks.map(([title, text]) =>
            el("div", { className: "mb-3 p-3 bg-white rounded shadow-sm" }, [
              el("h4", { className: "h6 fw-bold text-primary mb-1", text: title }),
              el("p", { className: "small text-muted mb-0", text })
            ])
          )
        ])
      ])
    ])
  ]);

  const children = [aboutLayout];
  analysis.figures.map(figureById).filter(Boolean).forEach((figure) => children.push(renderFigure(figure)));

  return el("div", { className: "d-flex flex-column gap-3" }, children);
}

function renderFigure(figure) {
  return el("div", { className: "card shadow-sm mb-4" }, [
    el("div", { className: "card-header bg-white d-flex justify-content-between align-items-start p-4 border-bottom-0" }, [
      el("div", {}, [
        el("span", { className: "badge bg-secondary mb-2", text: figure.section }),
        el("h3", { className: "h4 fw-bold mb-1", text: figure.title }),
        el("p", { className: "text-muted mb-0", text: figure.caption })
      ]),
      el("a", { className: "btn btn-outline-primary btn-sm", href: figure.src, target: "_blank", rel: "noreferrer", text: "Open full image" })
    ]),
    el("div", { className: "card-body text-center p-0 bg-light border-top" }, [
      el("img", { src: figure.src, alt: figure.title, className: "img-fluid p-3 figure-img-custom", loading: "lazy" })
    ])
  ]);
}

async function renderTable(table) {
  const parsed = await loadTableData(table);
  const search = state.search.trim().toLowerCase();
  const rows = parsed.records.filter((record) => {
    if (!search) return true;
    return Object.values(record).join(" ").toLowerCase().includes(search);
  });
  const shown = rows.slice(0, 220);
  
  const tableNode = el("table", { className: "table table-striped table-hover table-bordered mb-0" });
  const headRow = el("tr");
  parsed.columns.forEach((column) => headRow.appendChild(el("th", { text: column, className: "bg-light text-nowrap" })));
  tableNode.appendChild(el("thead", {}, [headRow]));
  
  const body = el("tbody");
  shown.forEach((record) => {
    const tr = el("tr");
    parsed.columns.forEach((column) => tr.appendChild(el("td", { text: compact(record[column]), className: "text-nowrap" })));
    body.appendChild(tr);
  });
  tableNode.appendChild(body);

  const panel = el("div", { className: "card shadow-sm mb-4" }, [
    el("div", { className: "card-header bg-white p-4" }, [
      el("div", { className: "d-flex justify-content-between align-items-start mb-3" }, [
        el("div", {}, [
          el("span", { className: "badge bg-secondary mb-2", text: table.section }),
          el("h3", { className: "h4 fw-bold mb-1", text: table.title }),
          el("p", { className: "text-muted mb-0", text: table.description })
        ]),
        el("a", { className: "btn btn-primary btn-sm", href: table.path, download: table.path.split("/").pop(), text: "Download CSV" })
      ]),
      el("div", { className: "row g-3 align-items-center" }, [
        el("div", { className: "col-md-6" }, [
          el("label", { className: "visually-hidden", text: "Search this table" }),
          el("input", { type: "search", className: "form-control", value: state.search, placeholder: "Search outcome, gene, locus..." })
        ]),
        el("div", { className: "col-md-6 text-md-end text-muted small" }, [
          el("span", { text: `${rows.length.toLocaleString()} matching rows; showing ${shown.length.toLocaleString()}.` })
        ])
      ])
    ]),
    el("div", { className: "card-body p-0 table-responsive" }, [tableNode])
  ]);

  panel.querySelector("input").addEventListener("input", (event) => {
    state.search = event.target.value;
    render();
  });
  return panel;
}

function renderDownloads() {
  const grid = el("div", { className: "row g-4 mt-3" });
  state.manifest.tables.forEach((table) => {
    grid.appendChild(el("div", { className: "col-md-6 col-lg-4" }, [
      el("div", { className: "card shadow-sm h-100" }, [
        el("div", { className: "card-body d-flex flex-column" }, [
          el("span", { className: "badge bg-secondary align-self-start mb-2", text: table.section }),
          el("h3", { className: "h5 fw-bold mb-2", text: table.title }),
          el("p", { className: "card-text text-muted small mb-4", text: table.description }),
          el("a", { className: "btn btn-outline-primary mt-auto", href: table.path, download: table.path.split("/").pop(), text: `Download ${table.format.toUpperCase()}` })
        ])
      ])
    ]));
  });
  return grid;
}

async function renderAnalysis(analysis) {
  if (analysis.id === "about") return renderAbout(analysis);
  
  if (analysis.id === "downloads") {
    return el("div", {}, [
      renderIntro(analysis),
      renderDownloads()
    ]);
  }

  const children = [renderIntro(analysis)];
  
  analysis.figures.map(figureById).filter(Boolean).forEach((figure) => children.push(renderFigure(figure)));
  
  for (const id of analysis.tables) {
    const table = tableById(id);
    if (table) children.push(await renderTable(table));
  }
  
  if (children.length === 1) {
    children.push(el("div", { className: "alert alert-secondary text-center mt-4", text: "No bundled outputs for this tab." }));
  }
  
  return el("div", { className: "d-flex flex-column gap-3" }, children);
}

function renderIntro(analysis) {
  return el("div", { className: "card shadow-sm mb-4 bg-primary text-white border-0" }, [
    el("div", { className: "card-body p-4" }, [
      el("span", { className: "badge bg-white text-primary mb-2", text: analysis.label }),
      el("h2", { className: "card-title fw-bold mb-2", text: analysis.title }),
      el("p", { className: "card-text lead mb-0 opacity-75", text: analysis.text })
    ])
  ]);
}

async function render() {
  renderTabs();
  const app = $("#app");
  const analysis = analyses.find((item) => item.id === state.activeTab) || analyses[0];
  app.innerHTML = "";
  app.appendChild(await renderAnalysis(analysis));
}

async function init() {
  const response = await fetch("manifest.json");
  state.manifest = await response.json();
  $("#updatedStamp").textContent = `Updated ${state.manifest.updated}`;
  await render();
}

init().catch((error) => {
  $("#app").innerHTML = "";
  $("#app").appendChild(el("div", { className: "alert alert-danger m-4" }, [
    el("h2", { className: "h4 alert-heading", text: "Browser failed to load" }),
    el("p", { className: "mb-0", text: error.message })
  ]));
});
