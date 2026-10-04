/** Tokens that mean "no value" in the proteomics software outputs this tool reads. */
export const MISSING_TOKENS = new Set(["", "nan", "na", "n/a", "#n/a", "filtered", "null", "n.d.", "nd", "-", "--", "inf", "-inf"]);

export const isMissingToken = (s: string | undefined) => s === undefined || MISSING_TOKENS.has(s.trim().toLowerCase());

/**
 * Parse one number. With decimalComma the comma is the decimal mark and points are thousands marks.
 * Otherwise a point is the decimal mark and commas are thousands marks. Returns NaN for anything else.
 */
export function parseNumber(s: string | undefined, decimalComma = false): number {
  if (isMissingToken(s)) return NaN;
  let t = s!.trim();
  if (decimalComma) t = t.replace(/\./g, "").replace(",", ".");
  else t = t.replace(/,/g, "");
  if (!/^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/.test(t)) return NaN;
  const v = Number(t);
  return Number.isFinite(v) ? v : NaN;
}

/**
 * Decide whether a set of cell values uses a decimal comma. Thousands separators show up as several comma
 * groups ("12,345,678") or as a point plus commas, which rules a decimal comma out. Values like "24,512" with a
 * single comma and no point are read as decimals, because raw intensities in a thousands style file almost always
 * have more than one comma group. Mixed or ambiguous input is treated as decimal point style.
 */
export function usesDecimalComma(values: (string | undefined)[]): boolean {
  let comma = 0, point = 0, thousands = 0;
  for (const v of values) {
    if (!v) continue;
    const t = v.trim();
    if (/^[-+]?\d{1,3}(,\d{3}){2,}(\.\d+)?$/.test(t) || /^[-+]?\d{1,3}(,\d{3})+\.\d+$/.test(t)) thousands++;
    else if (/^[-+]?\d+,\d+(e[-+]?\d+)?$/i.test(t)) comma++;
    else if (/\./.test(t)) point++;
  }
  return comma > 0 && point === 0 && thousands === 0;
}
