export interface Group {
  name: string;
  columns: string[];
}

export type Role = "bait" | "control" | "off";

const REPLICATE_SUFFIX = /[\s_.\-]*(?:rep(?:licate)?)?[\s_.\-]*\d+$/i;

/** Group LFQ columns that only differ by a trailing replicate number, keeping column order. */
export function groupColumns(columns: string[], clean: (column: string) => string = (c) => c.replace(/^LFQ intensity\s*/i, "")): Group[] {
  const groups = new Map<string, string[]>();
  for (const col of columns) {
    const stem = clean(col);
    const name = stem.replace(REPLICATE_SUFFIX, "") || stem;
    const list = groups.get(name);
    if (list) list.push(col); else groups.set(name, [col]);
  }
  return [...groups].map(([name, cols]) => ({ name, columns: cols }));
}

/** Suggest a role from the group name. Only obvious keywords; everything else stays off. */
export function suggestRole(name: string): Role {
  // An explicit "bait" wins over words like GFP or WT, so "GFP_bait" is a bait.
  if (/(^|[\s_.\-])bait([\s_.\-]|$)/i.test(name)) return "bait";
  if (/(^|[\s_.\-])(ctrl|control|neg|mock|igg|empty|gfp|wt|beads?)([\s_.\-]|$)/i.test(name)) return "control";
  return "off";
}
