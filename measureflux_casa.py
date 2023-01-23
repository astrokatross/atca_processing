# Script to be run WITHIN CASA!
day = "day4"
band="c"

filename = open("/home/cira/ATCA/bin/source_fluxesdict.json", "r")
srcs = json.load(filename)
src_names = list(srcs.keys())


for key in list(srcs.keys()):

    if band == "c":
        ms = f"c3487_{day}_cx.ms"
        spw = "0"
    elif band == "x":
        ms = f"c3487_{day}_cx.ms"
        spw = "1"
    elif band == "l":
        ms = f"c3487_{day}_l.ms"
        spw = ""
    try:
        dict_file = open(f"/home/cira/ATCA/bin/{key}_dict.json", "r")
        src_dict = json.load(dict_file)

    except:
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


        with open(f"/home/cira/ATCA/bin/{key}_dict.json", "w") as f: 
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

        dict_file = open(f"/home/cira/ATCA/bin/{fieldname}_dict.json", "r+")
        src_dict = json.load(dict_file)

        fitflux = src_dict[fieldname][f"{band}"][f"{day}"]

        uvmodelfit(
            vis=ms,
            niter=10,
            field=fieldname,
            selectdata=True,
            spw="0",
            scan = f"{scan}",
            sourcepar = [fitflux, 0, 0],
            outfile=f"{fieldname}_{band}_{day}_{scan}.cl"
        )

        tbl = cl.open(f"{fieldname}_{band}_{day}_{scan}.cl")
        fit = cl.getcomponent(0)
        flux = fit["flux"]["value"][0]
        src_dict[fieldname][f"{band}"][f"{timestamp}"] = flux

        with open(f"/home/cira/ATCA/bin/{fieldname}_dict.json", "w") as f: 
            json.dump(src_dict,f)
    else: 
        print("This is not a source I care about apparently!")