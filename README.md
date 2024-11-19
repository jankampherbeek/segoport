# segoport
A port of the Swiss Ephemeris to GoLang

This project is an attempt to port the code of the Swiss Ephemeris (SE) from C to Go. This would make it possible to use the SE without needing a FFI as cgo. That would benefit cross-platform use, especially for mobile devices, it would simplify accessing the SE and it would lead to a better performance.

The SE is large, about 60k lines of code. A comparable port for the C# language already exists: https://github.com/ygrenier/SwissEphNet. I will use the code of this project as a starting point. I will also use LLM's to assist in the conversion.

This is an ambitious project and there is no guarantee that I will succeed. But I will do my utmost.

The SeGoPort will be open source; the general GPL conditions apply. The original code from the SE is also covered by the GPL. If you want to write a commercial application, you can only do so if it is open source. You also need to take the license of the original SE into account. You can use the SE for closed source applications if you buy a commercial license. This port is however only available for open source applications.

Jan Kampherbeek

November, 11 2024.

