import { describe, expect, it } from "vitest";
import { groupColumns, suggestRole } from "../src/groups";

describe("groupColumns", () => {
  it("groups replicates by trailing number", () => {
    const g = groupColumns(["LFQ intensity BAIT_1", "LFQ intensity BAIT_2", "LFQ intensity CTRL_1", "LFQ intensity CTRL_2"]);
    expect(g.map((x) => x.name)).toEqual(["BAIT", "CTRL"]);
    expect(g[0].columns).toHaveLength(2);
  });
  it("handles messy separators and keeps order", () => {
    const g = groupColumns(["set01_Cxxc1___1", "set01_Cxxc1___2", "24-28hAPF_483_Ilk_01", "24-28hAPF_483_Ilk_02", "Sample rep3"]);
    expect(g.map((x) => x.name)).toEqual(["set01_Cxxc1", "24-28hAPF_483_Ilk", "Sample"]);
  });
  it("never produces an empty group name", () => {
    expect(groupColumns(["LFQ intensity 1", "LFQ intensity 2"]).every((x) => x.name.length > 0)).toBe(true);
  });
});

describe("suggestRole", () => {
  it("recognises obvious controls and baits only", () => {
    expect(suggestRole("CTRL")).toBe("control");
    expect(suggestRole("IgG_mock")).toBe("control");
    expect(suggestRole("BAIT")).toBe("bait");
    expect(suggestRole("Ilk")).toBe("off");
    expect(suggestRole("Controller")).toBe("off");
  });
});
