"""Query JASPAR REST API for every PPAR / RXR related profile (all collections, all versions)."""
import json, urllib.request, urllib.parse, time, sys, os
from concurrent.futures import ThreadPoolExecutor

BASE = "https://jaspar.elixir.no/api/v1"
OUT = os.path.dirname(os.path.abspath(__file__))


def get(url, tries=5):
    err = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json", "User-Agent": "Mozilla/5.0 (academic motif retrieval)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # noqa
            err = e
            time.sleep(2 + 2 * i)
    raise RuntimeError(f"failed {url}: {err}")


def paged(url):
    res = []
    while url:
        d = get(url)
        res.extend(d["results"])
        url = d.get("next")
    return res


collections = ["CORE", "CNEP", "FAM", "PBM", "PHYLOFACTS", "POLII", "SPLICE", "UNVALIDATED"]
terms = ["PPAR", "RXR"]

hits = {}
jobs = []
for t in terms:
    # default (no collection filter) plus each explicit collection
    jobs.append((t, None))
    for c in collections:
        jobs.append((t, c))


def run(job):
    t, c = job
    q = {"search": t, "page_size": 500, "format": "json"}
    if c:
        q["collection"] = c
    return job, paged(f"{BASE}/matrix/?" + urllib.parse.urlencode(q))


with ThreadPoolExecutor(max_workers=4) as ex:
    for job, res in ex.map(run, jobs):
        print(f"search={job[0]:5s} collection={str(job[1]):12s} -> {len(res)} profiles", flush=True)
        for r in res:
            hits[r["matrix_id"]] = r

# also ask for exact gene names, in case 'search' misses something
for nm in ["PPARA", "PPARD", "PPARG", "RXRA", "RXRB", "RXRG", "PPARA::RXRA", "PPARG::RXRA", "PPARD::RXRA"]:
    res = paged(f"{BASE}/matrix/?" + urllib.parse.urlencode({"name": nm, "page_size": 500, "format": "json"}))
    print(f"name={nm:12s} -> {len(res)} profiles", flush=True)
    for r in res:
        hits[r["matrix_id"]] = r

print("unique latest-listing profiles:", len(hits))
base_ids = sorted({h["matrix_id"].split(".")[0] for h in hits.values()})
print("unique base IDs:", len(base_ids))

# every version of every base id
all_ids = set(hits)


def versions(b):
    try:
        d = paged(f"{BASE}/matrix/{b}/versions/?format=json&page_size=200")
        return b, [x["matrix_id"] for x in d]
    except Exception as e:  # noqa
        return b, []


with ThreadPoolExecutor(max_workers=4) as ex:
    for b, ids in ex.map(versions, base_ids):
        all_ids.update(ids)

print("all versioned IDs:", len(all_ids))


def detail(mid):
    return mid, get(f"{BASE}/matrix/{mid}/?format=json")


details = {}
with ThreadPoolExecutor(max_workers=4) as ex:
    for mid, d in ex.map(detail, sorted(all_ids)):
        details[mid] = d

with open(os.path.join(OUT, "jaspar_details_raw.json"), "w") as fh:
    json.dump(details, fh, indent=1)
print("saved", len(details), "details")
