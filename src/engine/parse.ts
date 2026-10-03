export interface Table {
  columns: string[];
  rows: string[][];
}

/** Parse delimited text. Handles quoted fields and CRLF. */
export function parseDelimited(text: string, sep: string, quote = '"'): Table {
  const rows: string[][] = [];
  let row: string[] = [];
  let field = "";
  let inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (inQuotes) {
      if (ch === quote) {
        if (text[i + 1] === quote) { field += quote; i++; } else inQuotes = false;
      } else field += ch;
    } else if (quote && ch === quote && field === "") {
      inQuotes = true;
    } else if (ch === sep) {
      row.push(field); field = "";
    } else if (ch === "\n" || ch === "\r") {
      if (ch === "\r" && text[i + 1] === "\n") i++;
      row.push(field); field = "";
      if (row.length > 1 || row[0] !== "") rows.push(row);
      row = [];
    } else field += ch;
  }
  if (field !== "" || row.length) { row.push(field); rows.push(row); }
  const columns = rows.shift() ?? [];
  return { columns, rows };
}

export function sniffSeparator(text: string): string {
  const head = text.slice(0, 4096).split(/\r?\n/)[0] ?? "";
  return (head.match(/\t/g)?.length ?? 0) >= (head.match(/,/g)?.length ?? 0) ? "\t" : ",";
}
