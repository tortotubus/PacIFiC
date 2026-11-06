#include "Libxml2Backend.hh"

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

class TemporaryDirectory {
public:
    TemporaryDirectory()
    {
        auto pattern = (std::filesystem::temp_directory_path()
                        / "grains3dgpu-libxml2-XXXXXX")
                           .string();
        std::vector<char> writable(pattern.begin(), pattern.end());
        writable.push_back('\0');
        const char* directory = mkdtemp(writable.data());
        if(!directory)
            throw std::runtime_error("Cannot create XML test directory");
        path_ = directory;
    }

    ~TemporaryDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

TEST(Libxml2BackendTest, ExposesTypedOwningElementHandles)
{
    const auto document = grains::xml::Libxml2Backend{}.parseString(
        R"(<Root count="2"><Item id="3">4.5</Item><Item id="7">8.5</Item></Root>)");
    const auto root = document.root();

    EXPECT_EQ(root.name(), "Root");
    EXPECT_EQ(root.requireAttribute<int>("count"), 2);
    ASSERT_EQ(root.getChildren("Item").size(), 2u);
    EXPECT_EQ(root.getChildren("Item")[0].requireAttribute<unsigned>("id"), 3u);
    EXPECT_DOUBLE_EQ(root.getChildren("Item")[1].value<double>(), 8.5);
}

TEST(Libxml2BackendTest, LoadsDefaultsOnlyFromTheConfiguredDtdDirectory)
{
    TemporaryDirectory temporary;
    const auto dtd = temporary.path() / "defaults.dtd";
    {
        std::ofstream output(dtd);
        ASSERT_TRUE(output.good());
        output << "<!ELEMENT Root EMPTY>\n"
                  "<!ATTLIST Root mode CDATA \"from-dtd\">\n";
    }

    grains::xml::ParseOptions options;
    options.dtd = grains::xml::DtdPolicy::LoadDefaults;
    options.dtdDirectory = temporary.path();
    const auto document = grains::xml::Libxml2Backend{}.parseString(
        "<!DOCTYPE Root SYSTEM \"defaults.dtd\"><Root/>",
        temporary.path() / "input.xml",
        options);

    EXPECT_EQ(document.root().requireAttribute<std::string>("mode"), "from-dtd");
}

} // namespace
