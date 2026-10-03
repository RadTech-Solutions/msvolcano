declare module "plotly.js-cartesian-dist-min" {
  const Plotly: {
    react(el: HTMLElement, data: unknown[], layout?: unknown, config?: unknown): Promise<unknown>;
    purge(el: HTMLElement): void;
    downloadImage(el: HTMLElement, opts: { format: string; filename: string; width?: number; height?: number }): Promise<unknown>;
  };
  export default Plotly;
}
