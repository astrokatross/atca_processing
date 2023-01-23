# Script to be run WITHIN CASA!
day = "day4"
band="x"


filename = open("/home/cira/ATCA/bin/data/source_fluxesdict.json", "r")
srcs = json.load(filename)
src_names = list(srcs.keys())


if band == "c":
    ms = f"c3487_{day}_cx.ms"
    spw = "0"
elif band == "x":
    ms = f"c3487_{day}_cx.ms"
    spw = "1"
elif band == "l":
    ms = f"c3487_{day}_l.ms"
    spw = ""

# TODO: add option to check if the band is there not just the file since it crashes for second freq atm
for key in list(srcs.keys()):
    # try:
    #     dict_file = open(f"/home/cira/ATCA/bin/data/{key}_dict.json", "r")
    #     src_dict = json.load(dict_file)

    # except:
    print("Couldn't find the file, refitting the day uv")
    uvmodelfit(
        vis=ms,
        niter=10,
        field=key,
        selectdata=True,
        spw="0",
        outfile=f"{key}_{band}_{day}.cl"
    )
    tbl = cl.open(f"{key}_{band}_{day}.cl")
    fit = cl.getcomponent(0)
    flux = fit["flux"]["value"][0]
    srcs[key][f"{band}"][f"{day}"] = [flux]
    src_dict = {f"{key}": {f"{band}": {f"{day}": flux}}}
    # src_dict[f"{key}"][f"{band}"][f"{day}"] = [flux]


    with open(f"/home/cira/ATCA/bin/data/{key}_dict.json", "w") as f: 
        json.dump(src_dict,f)


obsinfo = listobs(f"{ms}")
info_keys = list(obsinfo.keys())
for key in info_keys:
    if key.startswith("scan"):
        continue
    else:
        obsinfo.pop(key)

for key in list(obsinfo.keys()):
    fieldname = obsinfo[key]["0"]["FieldName"]
    if fieldname in src_names:
        timestamp = obsinfo[key]["0"]["BeginTime"]
        scan = obsinfo[key]["0"]["scanId"]

        dict_file = open(f"/home/cira/ATCA/bin/data/{fieldname}_dict.json", "r+")
        src_dict = json.load(dict_file)

        fitflux = src_dict[fieldname][f"{band}"][f"{day}"]
        uvmodelfit(
            vis=ms,
            niter=10,
            field=fieldname,
            selectdata=True,
            spw=spw,
            scan = f"{scan}",
            sourcepar = [fitflux, 0, 0],
            outfile=f"{fieldname}_{band}_{day}_{scan}.cl"
        )

        tbl = cl.open(f"{fieldname}_{band}_{day}_{scan}.cl")
        fit = cl.getcomponent(0)
        flux = fit["flux"]["value"][0]
        src_dict[fieldname][f"{band}"][f"{timestamp}"] = flux

        with open(f"/home/cira/ATCA/bin/data/{fieldname}_dict.json", "w") as f: 
            json.dump(src_dict,f)
    else: 
        print("This is not a source I care about apparently!")