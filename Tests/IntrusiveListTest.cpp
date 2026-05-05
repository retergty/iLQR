#include <gtest/gtest.h>

#include "IntrusiveList.hpp"

namespace
{
struct TestNode : public IntrusiveListNode<TestNode>
{
    explicit TestNode(int value) : value(value) {}

    int value;
};
}

TEST(IntrusiveListTest, EmptyListBeginsAtEnd)
{
    IntrusiveList<TestNode> list;

    EXPECT_TRUE(list.begin());
    EXPECT_EQ(list.begin(), list.end());
}

TEST(IntrusiveListTest, InsertBeforeEndPreservesOrder)
{
    IntrusiveList<TestNode> list;
    TestNode first(1);
    TestNode second(2);
    TestNode third(3);

    EXPECT_EQ(list.insert(list.end(), first)->value, 1);
    list.insert(list.end(), second);
    list.insert(list.end(), third);

    auto it = list.begin();
    ASSERT_NE(it, list.end());
    EXPECT_EQ((*it).value, 1);
    EXPECT_EQ((it++)->value, 1);
    ASSERT_NE(it, list.end());
    EXPECT_EQ(it->value, 2);
    EXPECT_EQ((it + 1)->value, 3);

    it += 1;
    ASSERT_NE(it, list.end());
    EXPECT_EQ(it->value, 3);
    ++it;
    EXPECT_EQ(it, list.end());

    --it;
    EXPECT_EQ(it->value, 3);
    EXPECT_EQ((it--)->value, 3);
    EXPECT_EQ(it->value, 2);
}

TEST(IntrusiveListTest, InsertBeforeBeginUpdatesHead)
{
    IntrusiveList<TestNode> list;
    TestNode first(1);
    TestNode zeroth(0);

    list.insert(list.end(), first);
    list.insert(list.begin(), zeroth);

    auto it = list.begin();
    ASSERT_NE(it, list.end());
    EXPECT_EQ(it->value, 0);
    ++it;
    ASSERT_NE(it, list.end());
    EXPECT_EQ(it->value, 1);
    ++it;
    EXPECT_EQ(it, list.end());
}

TEST(IntrusiveListTest, ConstIteratorTraversesReadOnlyList)
{
    IntrusiveList<TestNode> list;
    TestNode first(1);
    TestNode second(2);

    list.insert(list.end(), first);
    list.insert(list.end(), second);

    const IntrusiveList<TestNode> &constList = list;
    auto it = constList.begin();
    ASSERT_NE(it, constList.end());
    EXPECT_EQ(it->value, 1);
    EXPECT_EQ((it++)->value, 1);
    ASSERT_NE(it, constList.end());
    EXPECT_EQ((*it).value, 2);
    EXPECT_EQ(constList.cbegin()->value, 1);
    EXPECT_EQ((constList.cbegin() + 1)->value, 2);
}

TEST(IntrusiveListTest, EraseReturnsFollowingIteratorAndAllowsReinsert)
{
    IntrusiveList<TestNode> list;
    IntrusiveList<TestNode> other;
    TestNode first(1);
    TestNode second(2);
    TestNode third(3);

    list.insert(list.end(), first);
    list.insert(list.end(), second);
    list.insert(list.end(), third);

    auto afterErased = list.erase(++list.begin());
    ASSERT_NE(afterErased, list.end());
    EXPECT_EQ(afterErased->value, 3);

    auto it = list.begin();
    ASSERT_NE(it, list.end());
    EXPECT_EQ(it->value, 1);
    ++it;
    ASSERT_NE(it, list.end());
    EXPECT_EQ(it->value, 3);
    ++it;
    EXPECT_EQ(it, list.end());

    other.insert(other.end(), second);
    ASSERT_NE(other.begin(), other.end());
    EXPECT_EQ(other.begin()->value, 2);
}

TEST(IntrusiveListTest, RemoveUpdatesOwnerHeadAndDetachedNodeCanMove)
{
    IntrusiveList<TestNode> source;
    IntrusiveList<TestNode> destination;
    TestNode first(1);
    TestNode second(2);
    TestNode third(3);

    source.insert(source.end(), first);
    source.insert(source.end(), second);
    destination.insert(destination.end(), third);

    destination.insert(destination.end(), first);

    ASSERT_NE(source.begin(), source.end());
    EXPECT_EQ(source.begin()->value, 2);
    EXPECT_EQ(++source.begin(), source.end());

    auto destIt = destination.begin();
    ASSERT_NE(destIt, destination.end());
    EXPECT_EQ(destIt->value, 3);
    ++destIt;
    ASSERT_NE(destIt, destination.end());
    EXPECT_EQ(destIt->value, 1);
    ++destIt;
    EXPECT_EQ(destIt, destination.end());

    EXPECT_EQ(source.remove(&second), 0);
    EXPECT_EQ(source.begin(), source.end());
    EXPECT_EQ(source.remove(&second), -1);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
