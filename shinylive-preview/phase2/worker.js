// Static resources only. Shinylive still handles app frames and Shiny traffic.
const pawehScope = self.registration.scope;
const pawehCachePrefix = "paweh-shinylive:" + pawehScope + ":";
const pawehCacheName = pawehCachePrefix + PAWEH_ASSET_REVISION;

function pawehAssetRequest(request) {
  if (request.method !== "GET" || request.headers.has("Range") ||
      !request.url.startsWith(pawehScope)) return false;
  return Object.hasOwn(PAWEH_ASSETS, request.url.slice(pawehScope.length));
}

async function pawehAsset(event) {
  const request = event.request;
  const relative = request.url.slice(pawehScope.length);
  const contentType = PAWEH_ASSETS[relative];
  let cache;
  try {
    cache = await caches.open(pawehCacheName);
  } catch {
    // A full or unavailable cache must not prevent the app from starting.
  }
  async function load(resource) {
    const cached = await cache?.match(resource);
    if (cached) return cached;
    const response = await fetch(resource);
    if (cache && response.ok) {
      event.waitUntil(cache.put(resource, response.clone()).catch(() => {}));
    }
    return response;
  }
  if (contentType && typeof DecompressionStream === "function") {
    try {
      const response = await load(new URL(relative + ".gz", pawehScope).href);
      if (response.ok && response.body) {
        // Some static hosts decode .gz at the HTTP layer. Inspect the actual
        // response bytes before applying a second layer of decompression.
        const reader = response.body.getReader();
        const first = await reader.read();
        const gzipped = first.value?.[0] === 31 && first.value?.[1] === 139;
        const body = new ReadableStream({
          start(controller) {
            if (first.value) controller.enqueue(first.value);
            if (first.done) controller.close();
          },
          async pull(controller) {
            const next = await reader.read();
            if (next.value) controller.enqueue(next.value);
            if (next.done) controller.close();
          },
          cancel(reason) { return reader.cancel(reason); }
        });
        return new Response(gzipped ? body.pipeThrough(new DecompressionStream("gzip")) : body, {
          headers: { "Content-Type": contentType }
        });
      }
    } catch {
      // Use the uncompressed asset if the encoded request is unavailable.
    }
  }
  return load(request);
}
