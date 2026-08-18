// TRANSCRIBED BY HAND from benchmark/results/summary.tsv (every value copied
// verbatim; feature counts are `wc -l` of benchmark/data/<ds>.features.tsv
// minus the header line).
//
// NOT generated: no script on disk writes this file, and the previous header
// said "GENERATED ... do not edit by hand", which discouraged the one check that
// would catch drift. The audit re-derived every cell and found them exact, so the
// content is right and the label was wrong. A generator is owed as part of the
// re-run's "re-derive every table by script" step (TODO.md), where it can be
// validated against fresh output rather than against numbers it was written from.
#table(
  columns: 8,
  align: (left, right, left, right, right, right, right, right),
  table.header([Dataset], [Features], [Engine], [1% FDR], [5% FDR],
    [user-CPU (s)], [wall (s)], [peak RSS (MiB)]),
  [09746], [59,645], [uno], [862], [1209], [350.8], [33.0], [3414],
  [], [], [ms1searchpy rep 1], [837], [1126], [1457.4], [205.8], [7096],
  [], [], [ms1searchpy rep 2], [853], [1122], [1405.2], [195.2], [7075],
  [], [], [ms1searchpy rep 3], [823], [1042], [1411.0], [195.2], [6976],
  table.hline(),
  [09766], [80,033], [uno], [875], [1100], [361.5], [33.1], [3415],
  [], [], [ms1searchpy rep 1], [913], [1062], [1671.7], [233.4], [7659],
  [], [], [ms1searchpy rep 2], [804], [977], [1651.0], [226.8], [7548],
  [], [], [ms1searchpy rep 3], [789], [989], [1675.7], [232.1], [7603],
  table.hline(),
  [09734], [71,581], [uno], [1035], [1245], [353.2], [31.9], [3317],
  [], [], [ms1searchpy rep 1], [770], [1073], [1436.9], [203.7], [7090],
  [], [], [ms1searchpy rep 2], [880], [1024], [1413.9], [195.4], [7184],
  [], [], [ms1searchpy rep 3], [951], [1181], [1434.0], [203.8], [7230],
  table.hline(),
  [09778], [86,778], [uno], [1222], [1653], [416.8], [40.2], [3980],
  [], [], [ms1searchpy rep 1], [1173], [1660], [1647.3], [240.1], [7610],
  [], [], [ms1searchpy rep 2], [1220], [1510], [1654.0], [235.2], [7732],
  [], [], [ms1searchpy rep 3], [1147], [1474], [1668.2], [239.1], [7670],
  table.hline(),
  [09950], [87,019], [uno], [1440], [1740], [375.8], [37.0], [3531],
  [], [], [ms1searchpy rep 1], [1169], [1645], [1678.5], [245.3], [7532],
  [], [], [ms1searchpy rep 2], [1469], [1751], [1640.1], [234.8], [7527],
  [], [], [ms1searchpy rep 3], [1279], [1740], [1676.5], [247.4], [7650],
  table.hline(),
  [LUM2_05433_VG], [100,949], [uno], [2072], [2916], [422.7], [44.6], [4200],
  [], [], [ms1searchpy rep 1], [2081], [2629], [1557.9], [243.4], [7905],
  [], [], [ms1searchpy rep 2], [1481], [2286], [1537.8], [234.4], [7877],
  [], [], [ms1searchpy rep 3], [1483], [2319], [1536.4], [232.9], [7584],
  table.hline(),
)
