3-Jul-15
3-Jul-15
21-Aug-15
21-Aug-15
13-Oct-15
13-Oct-15
7-Jan-16
7-Jan-16
12-Apr-16
12-Apr-16

ref = as.Date("21-Aug-15", format = "%d-%b-%y")
-63 - as.numeric(ref - as.Date("13-Oct-15", format = "%d-%b-%y"))

c1d1ip = as.Date("23-Jul-15", format = "%d-%b-%y")
as.numeric(as.Date("21-Aug-15", format = "%d-%b-%y") - c1d1ip) - 16*7
as.numeric(as.Date("3-Jul-15", format = "%d-%b-%y") - c1d1ip) - 16*7

as.numeric(as.Date("13-Oct-15", format = "%d-%b-%y") - c1d1ip) - 16*7
as.numeric(as.Date("7-Jan-16", format = "%d-%b-%y") - c1d1ip) - 16*7
as.numeric(as.Date("12-Apr-16", format = "%d-%b-%y") - c1d1ip) - 16*7


as.numeric(as.Date("13-Oct-15", format = "%d-%b-%y") - c1d1ip) - 16*7

zero = as.Date("13-Oct-15", format = "%d-%b-%y")
as.numeric(as.Date("21-Aug-15", format = "%d-%b-%y") - zero)
as.numeric(as.Date("3-Jul-15", format = "%d-%b-%y") - zero)
as.numeric(as.Date("7-Jan-16", format = "%d-%b-%y") - zero)
as.numeric(as.Date("12-Apr-16", format = "%d-%b-%y") - zero)


